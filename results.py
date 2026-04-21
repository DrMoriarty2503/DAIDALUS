#!/usr/bin/env python3
"""
Анализ зависимости вероятности невыявления конфликтов от уровня шума.
Строит графики для СКО и HFOM (2×СКО).
"""

import os
import sys
import re
import csv
import math
import random
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from scipy import stats
from scipy.optimize import curve_fit

# ============================================================
# 1. КОНСТАНТЫ
# ============================================================
LOWC_HORIZ_M = 609.9
LOWC_VERT_FT = 500.0
WARNING_WINDOW_BEFORE = 15

# Параметры Земли
EARTH_A = 6378137.0
EARTH_F = 1 / 298.257223563
EARTH_B = EARTH_A * (1 - EARTH_F)
EARTH_E2 = 2 * EARTH_F - EARTH_F ** 2
EARTH_B_A = math.sqrt(1 - EARTH_E2)

# Количество симуляций для каждого уровня шума
NUM_SIMULATIONS_PER_NOISE = 50  # Перебираем сценарии несколько раз со случайным шумом


# ============================================================
# 2. ГЕОДЕЗИЧЕСКИЕ ФУНКЦИИ
# ============================================================

def geodetic_to_reduced(lat_deg: float) -> float:
    φ_rad = math.radians(lat_deg)
    tan_β = EARTH_B_A * math.tan(φ_rad)
    return math.atan(tan_β)


def andoyer_distance(lat1_deg: float, lon1_deg: float,
                     lat2_deg: float, lon2_deg: float) -> float:
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
    """Читает .daa файл и возвращает траектории"""
    trajectories = {}

    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

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

    result = {}
    for t, data in trajectories.items():
        if 'Ownship' in data and 'AC1' in data:
            result[t] = {
                'ownship': data['Ownship'],
                'intruder': data['AC1']
            }

    return result


def has_warning_in_row(row: Dict) -> bool:
    """Проверяет наличие предупреждения в строке CSV"""
    alert_level = row.get('alert_level', 'None')
    if alert_level not in ['None', 'Unknown', '']:
        return True

    for key, value in row.items():
        if ('h_band_' in key or 'v_band_' in key) and '_type' in key:
            if value in ['WARNING', 'NEAR', 'RECOVERY']:
                return True

    return False


def read_recommendations(filepath: str) -> Dict[int, bool]:
    """Читает CSV файл с рекомендациями"""
    warnings = {}

    if not os.path.exists(filepath):
        return warnings

    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                time_val = int(float(row.get('time', 0)))
                warnings[time_val] = has_warning_in_row(row)
            except (ValueError, KeyError):
                continue

    return warnings


def add_noise_to_position(lat: float, lon: float, noise_std_m: float) -> Tuple[float, float]:
    """
    Добавляет гауссов шум к координатам.
    Шум задается в метрах, конвертируется в градусы.
    """
    # Конвертируем метры в градусы (приблизительно)
    lat_rad = math.radians(lat)
    std_lat_deg = noise_std_m / 111320.0
    std_lon_deg = noise_std_m / (111320.0 * math.cos(lat_rad))

    noisy_lat = lat + random.gauss(0, std_lat_deg)
    noisy_lon = lon + random.gauss(0, std_lon_deg)

    return noisy_lat, noisy_lon


def find_conflict_periods(trajectories: Dict[int, Dict]) -> List[Tuple[int, int, float, float]]:
    """Находит периоды конфликта"""
    conflict_periods = []
    in_conflict = False
    period_start = None
    min_horiz = float('inf')
    min_vert = float('inf')

    for t in sorted(trajectories.keys()):
        data = trajectories[t]
        ownship = data['ownship']
        intruder = data['intruder']

        horiz_dist = andoyer_distance(ownship[0], ownship[1], intruder[0], intruder[1])
        vert_dist = abs(ownship[2] - intruder[2])

        is_conflict = (horiz_dist < LOWC_HORIZ_M and vert_dist < LOWC_VERT_FT)

        if is_conflict and not in_conflict:
            in_conflict = True
            period_start = t
            min_horiz = horiz_dist
            min_vert = vert_dist
        elif is_conflict and in_conflict:
            min_horiz = min(min_horiz, horiz_dist)
            min_vert = min(min_vert, vert_dist)
        elif not is_conflict and in_conflict:
            conflict_periods.append((period_start, t - 1, min_horiz, min_vert))
            in_conflict = False
            min_horiz = float('inf')
            min_vert = float('inf')

    if in_conflict:
        conflict_periods.append((period_start, max(trajectories.keys()), min_horiz, min_vert))

    return conflict_periods


def check_detection_for_noise_level(daa_path: str, csv_path: str,
                                    noise_std_m: float) -> bool:
    """
    Проверяет, обнаружил ли DAIDALUS конфликт при данном уровне шума.
    """
    # Реальные траектории
    real_trajectories = read_daa_trajectory(daa_path)
    recommendations = read_recommendations(csv_path)

    if not real_trajectories or not recommendations:
        return True  # Нет данных - считаем как обнаружено

    # Реальные конфликты
    real_conflicts = find_conflict_periods(real_trajectories)

    if not real_conflicts:
        return True  # Нет реального конфликта

    # Симулируем зашумленные данные
    noisy_trajectories = {}
    for t, data in real_trajectories.items():
        ownship = data['ownship']
        intruder = data['intruder']

        noisy_ownship_lat, noisy_ownship_lon = add_noise_to_position(
            ownship[0], ownship[1], noise_std_m
        )
        noisy_intruder_lat, noisy_intruder_lon = add_noise_to_position(
            intruder[0], intruder[1], noise_std_m
        )

        noisy_trajectories[t] = {
            'ownship': (noisy_ownship_lat, noisy_ownship_lon, ownship[2]),
            'intruder': (noisy_intruder_lat, noisy_intruder_lon, intruder[2])
        }

    # Находим конфликты в зашумленных данных
    noisy_conflicts = find_conflict_periods(noisy_trajectories)

    # Проверяем предупреждение в окне
    first_conflict = real_conflicts[0]
    conflict_start = first_conflict[0]
    window_start = max(0, conflict_start - WARNING_WINDOW_BEFORE)

    has_warning = any(
        t in recommendations and recommendations[t]
        for t in range(window_start, conflict_start + 1)
    )

    # Обнаружение: есть конфликт в зашумленных данных И есть предупреждение
    detected = (len(noisy_conflicts) > 0) and has_warning

    return detected


# ============================================================
# 4. ОСНОВНОЙ АНАЛИЗ
# ============================================================

def analyze_noise_dependency(folder_path: str,
                             noise_levels: List[float],
                             num_simulations: int = 50) -> Dict:
    """
    Анализирует зависимость вероятности невыявления от уровня шума.
    """
    # Поиск всех пар файлов
    daa_files = {}
    csv_files = {}

    for f in os.listdir(folder_path):
        if f.endswith('.daa'):
            match = re.search(r'conflict_(\d+)', f)
            if match:
                daa_files[int(match.group(1))] = os.path.join(folder_path, f)
        elif f.endswith('.csv') and 'recommendations' in f:
            match = re.search(r'conflict_(\d+)_recommendations', f)
            if match:
                csv_files[int(match.group(1))] = os.path.join(folder_path, f)

    common_nums = sorted(set(daa_files.keys()) & set(csv_files.keys()))

    if not common_nums:
        print("❌ Не найдено пар файлов")
        return {}

    print(f"Найдено {len(common_nums)} сценариев с конфликтами")
    print(f"Анализируем {len(noise_levels)} уровней шума")
    print(f"По {num_simulations} симуляций на уровень\n")

    results = {}

    for noise_std in noise_levels:
        print(f"Обработка σ = {noise_std:.1f} м (HFOM = {2 * noise_std:.1f} м)...")

        undetected_count = 0
        total_detections = 0

        for scenario_num in common_nums:
            daa_path = daa_files[scenario_num]
            csv_path = csv_files[scenario_num]

            # Запускаем несколько симуляций для каждого сценария
            for sim in range(num_simulations):
                # Используем seed для воспроизводимости
                random.seed(scenario_num * 10000 + sim)

                detected = check_detection_for_noise_level(daa_path, csv_path, noise_std)
                total_detections += 1

                if not detected:
                    undetected_count += 1

        miss_rate = undetected_count / total_detections * 100
        results[noise_std] = {
            'miss_rate': miss_rate,
            'undetected': undetected_count,
            'total': total_detections
        }

        print(f"  → Невыявлено: {undetected_count}/{total_detections} ({miss_rate:.2f}%)")

    return results


# ============================================================
# 5. ПОСТРОЕНИЕ ГРАФИКОВ
# ============================================================

def plot_results(results: Dict, folder_path: str):
    """Строит графики зависимости"""

    noise_levels = sorted(results.keys())
    miss_rates = [results[n]['miss_rate'] for n in noise_levels]
    hfom_levels = [2 * n for n in noise_levels]

    # Настройка стиля
    plt.style.use('seaborn-v0_8-darkgrid')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Цвета и маркеры
    colors = plt.cm.viridis(np.linspace(0, 1, len(noise_levels)))

    # График 1: Зависимость от СКО (σ)
    ax1.plot(noise_levels, miss_rates, 'o-', linewidth=2, markersize=8,
             color='#1f77b4', label='Экспериментальные данные')

    # Аппроксимация сигмоидой (логистическая функция)
    def sigmoid(x, L, k, x0):
        return L / (1 + np.exp(-k * (x - x0)))

    try:
        x_data = np.array(noise_levels)
        y_data = np.array(miss_rates)

        # Нормализуем y_data для сигмоиды
        y_norm = y_data / 100.0

        popt, _ = curve_fit(sigmoid, x_data, y_norm,
                            p0=[1, 0.05, 500], maxfev=5000)

        x_smooth = np.linspace(min(noise_levels), max(noise_levels), 200)
        y_smooth = sigmoid(x_smooth, *popt) * 100

        ax1.plot(x_smooth, y_smooth, '--', linewidth=2, color='#ff7f0e',
                 label=f'Аппроксимация')

        # Находим точку 50% невыявления
        if min(y_norm) <= 0.5 <= max(y_norm):
            x_50 = np.interp(0.5, y_norm, x_data)
            ax1.axvline(x=x_50, color='red', linestyle=':', alpha=0.7)
            ax1.text(x_50 + 5, 55, f'σ_50% = {x_50:.1f} м',
                     fontsize=9, color='red')
    except:
        pass

    ax1.set_xlabel('СКО (σ), м', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Вероятность невыявления, %', fontsize=12, fontweight='bold')
    ax1.set_title('Зависимость невыявления LoWC от СКО', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best')
    ax1.set_ylim(-5, 105)

    # Добавляем точки с цветовой кодировкой
    for i, (noise, rate) in enumerate(zip(noise_levels, miss_rates)):
        ax1.annotate(f'{rate:.1f}%', (noise, rate),
                     textcoords="offset points", xytext=(0, 10),
                     ha='center', fontsize=8)

    # График 2: Зависимость от HFOM (2×σ)
    ax2.plot(hfom_levels, miss_rates, 's-', linewidth=2, markersize=8,
             color='#2ca02c', label='Экспериментальные данные')

    # Аппроксимация для HFOM
    try:
        popt_hfom, _ = curve_fit(sigmoid, np.array(hfom_levels), y_norm,
                                 p0=[1, 0.025, 1000], maxfev=5000)

        x_hfom_smooth = np.linspace(min(hfom_levels), max(hfom_levels), 200)
        y_hfom_smooth = sigmoid(x_hfom_smooth, *popt_hfom) * 100

        ax2.plot(x_hfom_smooth, y_hfom_smooth, '--', linewidth=2, color='#d62728',
                 label=f'Аппроксимация')

        # Находим точку 50% невыявления
        if min(y_norm) <= 0.5 <= max(y_norm):
            x_hfom_50 = np.interp(0.5, y_norm, hfom_levels)
            ax2.axvline(x=x_hfom_50, color='red', linestyle=':', alpha=0.7)
            ax2.text(x_hfom_50 + 20, 55, f'HFOM_50% = {x_hfom_50:.1f} м',
                     fontsize=9, color='red')
    except:
        pass

    ax2.set_xlabel('HFOM (2×σ), м', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Вероятность невыявления, %', fontsize=12, fontweight='bold')
    ax2.set_title('Зависимость невыявления LoWC от HFOM', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')
    ax2.set_ylim(-5, 105)

    # Добавляем точки
    for i, (hfom, rate) in enumerate(zip(hfom_levels, miss_rates)):
        ax2.annotate(f'{rate:.1f}%', (hfom, rate),
                     textcoords="offset points", xytext=(0, 10),
                     ha='center', fontsize=8)

    plt.tight_layout()

    # Сохраняем графики
    plot_path = os.path.join(folder_path, 'noise_dependency_analysis.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"\n📊 График сохранен: {plot_path}")

    plt.show()

    return fig


def save_results_table(results: Dict, folder_path: str):
    """Сохраняет таблицу результатов"""

    noise_levels = sorted(results.keys())

    table_path = os.path.join(folder_path, 'noise_analysis_results.csv')
    with open(table_path, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['СКО (σ), м', 'HFOM (2σ), м', 'Вероятность невыявления, %',
                         'Невыявлено', 'Всего'])

        for noise in noise_levels:
            data = results[noise]
            writer.writerow([
                f"{noise:.1f}",
                f"{2 * noise:.1f}",
                f"{data['miss_rate']:.2f}",
                data['undetected'],
                data['total']
            ])

    print(f"📄 Таблица результатов сохранена: {table_path}")

    # Вывод в консоль
    print("\n" + "=" * 70)
    print("РЕЗУЛЬТАТЫ АНАЛИЗА")
    print("=" * 70)
    print(f"{'СКО (м)':<12} {'HFOM (м)':<12} {'Невыявление (%)':<20} {'Статус':<15}")
    print("-" * 70)

    for noise in noise_levels:
        data = results[noise]
        hfom = 2 * noise
        miss_rate = data['miss_rate']

        if miss_rate < 10:
            status = "✅ Отлично"
        elif miss_rate < 30:
            status = "⚠️ Хорошо"
        elif miss_rate < 50:
            status = "⚠️ Средне"
        else:
            status = "❌ Плохо"

        print(f"{noise:<12.1f} {hfom:<12.1f} {miss_rate:<20.2f} {status:<15}")

    print("=" * 70)


# ============================================================
# 6. ТОЧКА ВХОДА
# ============================================================

def main():
    if len(sys.argv) != 2:
        print("Использование: python noise_analysis.py <папка_с_файлами>")
        print("Пример: python noise_analysis.py random_conflicts")
        sys.exit(1)

    folder_path = sys.argv[1]

    if not os.path.isdir(folder_path):
        print(f"Ошибка: папка не найдена — {folder_path}")
        sys.exit(1)

    # Уровни шума для анализа (в метрах)
    # Типичные значения: от 0 до 2000 м с шагом 50-100 м
    noise_levels = list(range(0, 501, 50)) + list(range(600, 2001, 100))

    # Количество симуляций для каждого уровня шума
    num_simulations = NUM_SIMULATIONS_PER_NOISE

    print("\n" + "=" * 70)
    print(" АНАЛИЗ ЗАВИСИМОСТИ НЕВЫЯВЛЕНИЯ LoWC ОТ УРОВНЯ ШУМА")
    print("=" * 70)
    print(f"Порог LoWC: {LOWC_HORIZ_M} м / {LOWC_VERT_FT} фт")
    print(f"Окно предупреждения: {WARNING_WINDOW_BEFORE} сек")
    print(f"Анализируемые СКО: от {min(noise_levels)} до {max(noise_levels)} м")
    print(f"Симуляций на уровень: {num_simulations}")
    print("=" * 70 + "\n")

    # Запуск анализа
    results = analyze_noise_dependency(folder_path, noise_levels, num_simulations)

    if not results:
        print("❌ Анализ не выполнен")
        sys.exit(1)

    # Сохранение результатов
    save_results_table(results, folder_path)

    # Построение графиков
    plot_results(results, folder_path)

    # Дополнительная статистика
    print("\n📈 ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ:")
    print("-" * 50)
    print("• HFOM = 2×СКО — стандартный показатель точности GNSS")
    print("• Рекомендуемые требования к HFOM для DAA:")
    print("  - HFOM < 100 м: Отличная точность")
    print("  - HFOM 100-300 м: Хорошая точность")
    print("  - HFOM 300-600 м: Удовлетворительная точность")
    print("  - HFOM > 600 м: Низкая точность, риск невыявления")


if __name__ == "__main__":
    main()