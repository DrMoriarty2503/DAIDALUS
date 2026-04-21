#!/usr/bin/env python3
"""
Валидатор выявления конфликтов LoWC для DAIDALUS.
Сравнивает файлы сценариев (.daa) и файлы рекомендаций (.csv) для оценки
эффективности системы предотвращения столкновений.
"""

import os
import sys
import re
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# ============================================================
# 1. КОНСТАНТЫ
# ============================================================
# Пороги LoWC (Loss of Well Clear)
LOWC_HORIZ_M = 609.9  # Горизонтальный порог в метрах
LOWC_VERT_FT = 500.0  # Вертикальный порог в футах

# Окно предупреждения: за сколько секунд до конфликта должно быть предупреждение
WARNING_WINDOW_BEFORE = 15  # секунд

# Параметры Земли для геодезических расчетов (WGS-84)
EARTH_A = 6378137.0
EARTH_F = 1 / 298.257223563
EARTH_B = EARTH_A * (1 - EARTH_F)
EARTH_E2 = 2 * EARTH_F - EARTH_F ** 2
EARTH_B_A = math.sqrt(1 - EARTH_E2)


# ============================================================
# 2. ГЕОДЕЗИЧЕСКИЕ ФУНКЦИИ (из генератора сценариев)
# ============================================================

def geodetic_to_reduced(lat_deg: float) -> float:
    """Геодезическая → приведенная широта"""
    φ_rad = math.radians(lat_deg)
    tan_β = EARTH_B_A * math.tan(φ_rad)
    return math.atan(tan_β)


def andoyer_distance(lat1_deg: float, lon1_deg: float,
                     lat2_deg: float, lon2_deg: float) -> float:
    """
    Расстояние между точками по поверхности эллипсоида (метры)
    Использует метод Андуайе для высокой точности
    """
    u1 = geodetic_to_reduced(lat1_deg)
    u2 = geodetic_to_reduced(lat2_deg)
    delta_lambda = math.radians(lon2_deg - lon1_deg)

    sin_u1, sin_u2 = math.sin(u1), math.sin(u2)
    cos_u1, cos_u2 = math.cos(u1), math.cos(u2)
    sin_dl, cos_dl = math.sin(delta_lambda), math.cos(delta_lambda)

    p = cos_u2 * sin_dl
    q = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_dl
    n = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_dl

    sigma = math.atan2(math.sqrt(p ** 2 + q ** 2), n)

    if sigma < 1e-12:
        return 0.0

    sin_sigma, cos_sigma = math.sin(sigma), math.cos(sigma)

    if abs(1 + cos_sigma) < 1e-12:
        M = 0
    else:
        M = (sigma - 3 * sin_sigma) / (1 + cos_sigma)

    if abs(1 - cos_sigma) < 1e-12:
        N = 0
    else:
        N = (sigma + 3 * sin_sigma) / (1 - cos_sigma)

    U = (sin_u1 + sin_u2) ** 2
    V = (sin_u1 - sin_u2) ** 2

    delta_sigma = -0.25 * EARTH_E2 * (M * U + N * V)
    return EARTH_A * (sigma + delta_sigma)


# ============================================================
# 3. ЧТЕНИЕ ФАЙЛОВ
# ============================================================

def read_daa_trajectory(filepath: str) -> Dict[int, Dict]:
    """
    Читает .daa файл и возвращает траектории обоих ВС.

    Формат файла:
    NAME, lat, lon, alt, vx, vy, vz, time
    unitless, [deg], [deg], [ft], [knot], [knot], [fpm], [s]

    Возвращает:
        dict: {time: {'ownship': (lat, lon, alt), 'intruder': (lat, lon, alt)}}
    """
    trajectories = {}

    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Пропускаем заголовки (первые 2 строки)
    for line in lines[2:]:
        line = line.strip()
        if not line:
            continue

        parts = [p.strip() for p in line.split(',')]
        if len(parts) != 8:
            continue

        name = parts[0]
        try:
            lat = float(parts[1])
            lon = float(parts[2])
            alt = float(parts[3])
            time_val = int(float(parts[7]))
        except ValueError:
            continue

        if time_val not in trajectories:
            trajectories[time_val] = {}

        trajectories[time_val][name] = (lat, lon, alt)

    # Преобразуем в удобный формат
    result = {}
    for t, data in trajectories.items():
        if 'Ownship' in data and 'AC1' in data:
            result[t] = {
                'ownship': data['Ownship'],
                'intruder': data['AC1']
            }

    return result


def has_warning_in_row(row: Dict) -> bool:
    """
    Проверяет, содержит ли строка рекомендаций предупреждение.

    Предупреждение считается активным если:
    1. alert_level не 'None' и не 'Unknown'
    2. ИЛИ есть горизонтальные или вертикальные bands типа 'WARNING', 'NEAR', 'RECOVERY'
    """
    # Проверяем alert_level
    alert_level = row.get('alert_level', 'None')
    if alert_level not in ['None', 'Unknown', '']:
        return True

    # Проверяем горизонтальные bands
    for key, value in row.items():
        if ('h_band_' in key or 'v_band_' in key) and '_type' in key:
            if value in ['WARNING', 'NEAR', 'RECOVERY']:
                return True

    return False


def read_recommendations(filepath: str) -> Dict[int, bool]:
    """
    Читает CSV файл с рекомендациями от DAIDALUS.

    Возвращает:
        dict: {time: has_warning}
    """
    warnings = {}

    if not os.path.exists(filepath):
        print(f"   ⚠️ Файл не найден: {filepath}")
        return warnings

    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            try:
                time_val = int(float(row.get('time', 0)))
                has_warning = has_warning_in_row(row)
                warnings[time_val] = has_warning

            except (ValueError, KeyError) as e:
                continue

    return warnings


# ============================================================
# 4. АНАЛИЗ КОНФЛИКТОВ
# ============================================================

def find_conflict_periods(trajectories: Dict[int, Dict]) -> List[Tuple[int, int, float, float]]:
    """
    Находит все периоды времени, когда ВС находятся в состоянии LoWC.

    LoWC = горизонтальное расстояние < 609.9м И вертикальное расстояние < 500фт

    Возвращает:
        список кортежей (start_time, end_time, min_horiz, min_vert)
    """
    conflict_periods = []
    in_conflict = False
    period_start = None
    min_horiz = float('inf')
    min_vert = float('inf')

    for t in sorted(trajectories.keys()):
        data = trajectories[t]
        ownship = data['ownship']
        intruder = data['intruder']

        # Рассчитываем расстояния
        horiz_dist = andoyer_distance(ownship[0], ownship[1], intruder[0], intruder[1])
        vert_dist = abs(ownship[2] - intruder[2])

        is_conflict = (horiz_dist < LOWC_HORIZ_M and vert_dist < LOWC_VERT_FT)

        if is_conflict and not in_conflict:
            # Начало конфликта
            in_conflict = True
            period_start = t
            min_horiz = horiz_dist
            min_vert = vert_dist
        elif is_conflict and in_conflict:
            # Продолжение конфликта - обновляем минимумы
            min_horiz = min(min_horiz, horiz_dist)
            min_vert = min(min_vert, vert_dist)
        elif not is_conflict and in_conflict:
            # Конец конфликта
            conflict_periods.append((period_start, t - 1, min_horiz, min_vert))
            in_conflict = False
            min_horiz = float('inf')
            min_vert = float('inf')

    # Если конфликт до конца симуляции
    if in_conflict:
        conflict_periods.append((period_start, max(trajectories.keys()), min_horiz, min_vert))

    return conflict_periods


def check_warning_before_conflict(recommendations: Dict[int, bool],
                                  conflict_start: int) -> bool:
    """
    Проверяет, было ли предупреждение за WARNING_WINDOW_BEFORE секунд до конфликта.

    Возвращает:
        True если предупреждение было, False если нет
    """
    window_start = max(0, conflict_start - WARNING_WINDOW_BEFORE)

    # Проверяем все моменты времени в окне
    for t in range(window_start, conflict_start + 1):
        if t in recommendations and recommendations[t]:
            return True

    return False


def analyze_single_pair(daa_path: str, csv_path: str) -> Dict:
    """
    Анализирует одну пару файлов (сценарий + рекомендации).

    Возвращает:
        dict: результаты анализа
    """
    # Читаем файлы
    trajectories = read_daa_trajectory(daa_path)
    recommendations = read_recommendations(csv_path)

    if not trajectories:
        return {
            'scenario': Path(daa_path).stem,
            'error': 'Нет данных траектории',
            'has_conflict': False,
            'detected': False
        }

    if not recommendations:
        return {
            'scenario': Path(daa_path).stem,
            'error': 'Нет данных рекомендаций',
            'has_conflict': False,
            'detected': False
        }

    # Находим периоды конфликтов
    conflict_periods = find_conflict_periods(trajectories)

    if not conflict_periods:
        return {
            'scenario': Path(daa_path).stem,
            'has_conflict': False,
            'detected': True,  # Нет конфликта - нечего выявлять
            'conflict_start': None,
            'conflict_end': None,
            'min_horiz_m': None,
            'min_vert_ft': None,
            'num_conflicts': 0,
            'warning_before': None,
            'error': None
        }

    # Проверяем каждый конфликтный период
    all_detected = True
    undetected_periods = []

    for start_t, end_t, min_h, min_v in conflict_periods:
        has_warning = check_warning_before_conflict(recommendations, start_t)
        if not has_warning:
            all_detected = False
            undetected_periods.append((start_t, end_t))

    # Берем первый конфликт для деталей
    first_conflict = conflict_periods[0]

    return {
        'scenario': Path(daa_path).stem,
        'has_conflict': True,
        'detected': all_detected,
        'conflict_start': first_conflict[0],
        'conflict_end': first_conflict[1],
        'min_horiz_m': first_conflict[2],
        'min_vert_ft': first_conflict[3],
        'num_conflicts': len(conflict_periods),
        'undetected_periods': undetected_periods,
        'error': None
    }


# ============================================================
# 5. АНАЛИЗ ПАПКИ
# ============================================================

def analyze_folder(folder_path: str) -> Dict:
    """
    Анализирует все пары файлов в папке.

    Ожидает структуру:
    - conflict_XXXX.daa - файлы сценариев
    - conflict_XXXX_recommendations.csv - файлы рекомендаций
    """

    print("\n" + "=" * 80)
    print(" АНАЛИЗ ЭФФЕКТИВНОСТИ ОБНАРУЖЕНИЯ LoWC")
    print("=" * 80)
    print(f"Пороги LoWC: {LOWC_HORIZ_M}м / {LOWC_VERT_FT}фт")
    print(f"Окно предупреждения: {WARNING_WINDOW_BEFORE} секунд ДО конфликта\n")

    # Поиск файлов
    daa_files = {}
    csv_files = {}

    for f in os.listdir(folder_path):
        filepath = os.path.join(folder_path, f)

        if f.endswith('.daa'):
            # Ищем номер сценария conflict_XXXX.daa
            match = re.search(r'conflict_(\d+)', f)
            if match:
                num = int(match.group(1))
                daa_files[num] = filepath

        elif f.endswith('.csv') and 'recommendations' in f:
            # Ищем номер сценария conflict_XXXX_recommendations.csv
            match = re.search(r'conflict_(\d+)_recommendations', f)
            if match:
                num = int(match.group(1))
                csv_files[num] = filepath

    # Находим общие сценарии
    common_nums = set(daa_files.keys()) & set(csv_files.keys())

    if not common_nums:
        print("❌ Не найдено пар .daa + _recommendations.csv с совпадающими номерами")
        print(f"   Найдено .daa: {sorted(daa_files.keys())}")
        print(f"   Найдено .csv: {sorted(csv_files.keys())}")
        return {}

    print(f"✅ Найдено {len(common_nums)} пар файлов для анализа\n")

    # Анализ каждого сценария
    results = []
    total_conflicts = 0
    detected_conflicts = 0
    undetected_conflicts = 0

    for num in sorted(common_nums):
        daa_path = daa_files[num]
        csv_path = csv_files[num]

        result = analyze_single_pair(daa_path, csv_path)
        results.append(result)

        # Вывод результатов
        if result.get('has_conflict'):
            total_conflicts += 1
            if result.get('detected'):
                detected_conflicts += 1
                print(f"✅ Сценарий {num:04d}: КОНФЛИКТ {result['conflict_start']}-{result['conflict_end']}с, "
                      f"предупреждение БЫЛО за {WARNING_WINDOW_BEFORE}с")
                print(f"      Мин. расстояния: H={result['min_horiz_m']:.1f}м, V={result['min_vert_ft']:.1f}фт")
            else:
                undetected_conflicts += 1
                print(f"❌ Сценарий {num:04d}: КОНФЛИКТ {result['conflict_start']}-{result['conflict_end']}с, "
                      f"НЕТ ПРЕДУПРЕЖДЕНИЯ за {WARNING_WINDOW_BEFORE}с")
                print(f"      Мин. расстояния: H={result['min_horiz_m']:.1f}м, V={result['min_vert_ft']:.1f}фт")
                if result.get('undetected_periods'):
                    print(f"      Необнаруженные периоды: {result['undetected_periods']}")
        else:
            if result.get('error'):
                print(f"⚠️  Сценарий {num:04d}: {result['error']}")
            else:
                print(f"ℹ️  Сценарий {num:04d}: Конфликтов LoWC не обнаружено")

    # Статистика
    detection_rate = (detected_conflicts / total_conflicts * 100) if total_conflicts > 0 else 0
    undetected_rate = (undetected_conflicts / total_conflicts * 100) if total_conflicts > 0 else 0

    stats = {
        'total_scenarios': len(results),
        'conflicts_found': total_conflicts,
        'detected_conflicts': detected_conflicts,
        'undetected_conflicts': undetected_conflicts,
        'detection_rate_percent': detection_rate,
        'undetected_rate_percent': undetected_rate,
        'warning_window_before': WARNING_WINDOW_BEFORE,
        'lowc_thresholds': {'horizontal_m': LOWC_HORIZ_M, 'vertical_ft': LOWC_VERT_FT}
    }

    # Вывод статистики
    print("\n" + "=" * 80)
    print(" СТАТИСТИКА ВАЛИДАЦИИ")
    print("=" * 80)
    print(f"Всего проанализировано сценариев:     {stats['total_scenarios']}")
    print(f"Сценариев с конфликтом LoWC:          {stats['conflicts_found']}")
    print(f"✅ Выявленных конфликтов:              {stats['detected_conflicts']}")
    print(f"❌ НЕвыявленных конфликтов:            {stats['undetected_conflicts']}")
    print(f"\n📊 Процент выявления (Detection Rate):  {stats['detection_rate_percent']:.2f}%")
    print(f"📊 Процент невыявления (Miss Rate):     {stats['undetected_rate_percent']:.2f}%")
    print(f"\nПараметры:")
    print(f"  - Порог LoWC (горизонт): {LOWC_HORIZ_M} м")
    print(f"  - Порог LoWC (вертикаль): {LOWC_VERT_FT} фт")
    print(f"  - Окно предупреждения:    {WARNING_WINDOW_BEFORE} секунд ДО конфликта")

    # Сохраняем детальный отчет
    report_path = os.path.join(folder_path, "lowc_validation_report.csv")
    with open(report_path, 'w', newline='', encoding='utf-8-sig') as f:
        fieldnames = ['scenario', 'has_conflict', 'detected', 'conflict_start', 'conflict_end',
                      'min_horiz_m', 'min_vert_ft', 'num_conflicts', 'error']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            writer.writerow({
                'scenario': r.get('scenario', ''),
                'has_conflict': r.get('has_conflict', False),
                'detected': r.get('detected', False),
                'conflict_start': r.get('conflict_start', ''),
                'conflict_end': r.get('conflict_end', ''),
                'min_horiz_m': f"{r['min_horiz_m']:.2f}" if r.get('min_horiz_m') else '',
                'min_vert_ft': f"{r['min_vert_ft']:.1f}" if r.get('min_vert_ft') else '',
                'num_conflicts': r.get('num_conflicts', ''),
                'error': r.get('error', '')
            })

    print(f"\n📄 Детальный отчет сохранен: {report_path}")

    # Сохраняем краткий отчет
    summary_path = os.path.join(folder_path, "lowc_validation_summary.txt")
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("ОТЧЕТ О ВАЛИДАЦИИ DAIDALUS (LoWC)\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Всего проанализировано сценариев:     {stats['total_scenarios']}\n")
        f.write(f"Сценариев с конфликтом LoWC:          {stats['conflicts_found']}\n")
        f.write(f"Выявленных конфликтов:                {stats['detected_conflicts']}\n")
        f.write(f"НЕвыявленных конфликтов:              {stats['undetected_conflicts']}\n")
        f.write(f"\nПроцент выявления (Detection Rate):  {stats['detection_rate_percent']:.2f}%\n")
        f.write(f"Процент невыявления (Miss Rate):     {stats['undetected_rate_percent']:.2f}%\n\n")
        f.write("Параметры:\n")
        f.write(f"  - Порог LoWC (горизонт): {LOWC_HORIZ_M} м\n")
        f.write(f"  - Порог LoWC (вертикаль): {LOWC_VERT_FT} фт\n")
        f.write(f"  - Окно предупреждения:    {WARNING_WINDOW_BEFORE} секунд ДО конфликта\n")

    print(f"📄 Краткий отчет сохранен: {summary_path}")

    return stats


# ============================================================
# 6. ТОЧКА ВХОДА
# ============================================================

def main():
    if len(sys.argv) != 2:
        print("Использование: python validator.py <папка_с_файлами>")
        print("Пример: python validator.py random_conflicts")
        print("\nОжидаемая структура файлов:")
        print("  - conflict_0001.daa")
        print("  - conflict_0001_recommendations.csv")
        print("  - conflict_0002.daa")
        print("  - conflict_0002_recommendations.csv")
        print("  - ...")
        sys.exit(1)

    folder_path = sys.argv[1]

    if not os.path.isdir(folder_path):
        print(f"Ошибка: папка не найдена — {folder_path}")
        sys.exit(1)

    stats = analyze_folder(folder_path)

    if not stats:
        sys.exit(1)

    # Возвращаем код завершения
    # 0 - все конфликты выявлены, 1 - есть невыявленные
    if stats.get('undetected_conflicts', 0) > 0:
        print(f"\n❌ ВАЛИДАЦИЯ НЕ ПРОЙДЕНА: {stats['undetected_conflicts']} конфликтов не выявлено")
        sys.exit(1)
    else:
        print(f"\n✅ ВАЛИДАЦИЯ ПРОЙДЕНА: Все {stats['conflicts_found']} конфликтов выявлены")
        sys.exit(0)


if __name__ == "__main__":
    main()