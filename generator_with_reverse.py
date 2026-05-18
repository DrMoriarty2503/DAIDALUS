import math
import random
import os
import json
from collections import deque

from models import Ownship, Intruder
from daa_logic_docker import *

'''
    Генератор конфликтных сценариев с плавным уклонением и возвратом на линию пути
'''

# ============================================================
# 1. ПАРАМЕТРЫ WGS-84
# ============================================================
EARTH_A = 6378137.0
EARTH_F = 1 / 298.257223563
EARTH_B = EARTH_A * (1 - EARTH_F)
EARTH_E2 = 2 * EARTH_F - EARTH_F ** 2
EARTH_B_A = math.sqrt(1 - EARTH_E2)

type_conflict = "NMAC"
if type_conflict == "LOWC":
    SAFE_HORIZ_M = 609.9
    SAFE_VERT_FT = 500.0
else:
    SAFE_HORIZ_M = 150.0
    SAFE_VERT_FT = 100.0

KNOT_TO_MS = 0.514444

# Параметры плавности
HEADING_SMOOTHING = 0.08  # Плавность изменения курса при уклонении
RETURN_SMOOTHING = 0.03  # Плавность возврата к исходному курсу
RETURN_DISTANCE_M = 8000  # Расстояние до точки возврата


# ============================================================
# 2. ГЕОДЕЗИЧЕСКИЕ ФУНКЦИИ
# ============================================================

def heading_to_vx_vy(heading_deg, speed_knot):
    rad = math.radians(heading_deg)
    return speed_knot * math.sin(rad), speed_knot * math.cos(rad)


def vx_vy_to_heading(vx, vy):
    gs = math.sqrt(vx ** 2 + vy ** 2)
    heading_deg = math.degrees(math.atan2(vx, vy))
    return heading_deg if heading_deg >= 0 else heading_deg + 360, gs


def geodetic_to_reduced(lat_deg):
    return math.atan(EARTH_B_A * math.tan(math.radians(lat_deg)))


def andoyer_distance(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
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


def andoyer_azimuth(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    B1, B2 = math.radians(lat1_deg), math.radians(lat2_deg)
    delta_L = math.radians(lon2_deg - lon1_deg)
    y = math.cos(B2) * math.sin(delta_L)
    x = math.cos(B1) * math.sin(B2) - math.sin(B1) * math.cos(B2) * math.cos(delta_L)
    azimuth_deg = math.degrees(math.atan2(y, x))
    return azimuth_deg if azimuth_deg >= 0 else azimuth_deg + 360


def direct_geodetic_problem(lat1_deg, lon1_deg, azimuth_deg, distance_m):
    lat1 = math.radians(lat1_deg)
    lon1 = math.radians(lon1_deg)
    alpha1 = math.radians(azimuth_deg)

    U1 = math.atan((1 - EARTH_F) * math.tan(lat1))
    sin_U1, cos_U1 = math.sin(U1), math.cos(U1)
    sigma1 = math.atan2(math.tan(U1), math.cos(alpha1))
    sin_alpha = cos_U1 * math.sin(alpha1)
    cos2_alpha = 1 - sin_alpha ** 2
    u2 = cos2_alpha * (EARTH_A ** 2 - EARTH_B ** 2) / (EARTH_B ** 2)

    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    sigma = distance_m / (EARTH_B * A)

    for _ in range(100):
        cos_2sigma_m = math.cos(2 * sigma1 + sigma)
        sin_sigma, cos_sigma = math.sin(sigma), math.cos(sigma)
        delta_sigma = B * sin_sigma * (cos_2sigma_m + B / 4 * (cos_sigma * (-1 + 2 * cos_2sigma_m ** 2) -
                                                               B / 6 * cos_2sigma_m * (-3 + 4 * sin_sigma ** 2) * (
                                                                       -3 + 4 * cos_2sigma_m ** 2)))
        sigma_new = distance_m / (EARTH_B * A) + delta_sigma
        if abs(sigma_new - sigma) < 1e-12:
            sigma = sigma_new
            break
        sigma = sigma_new

    sin_sigma, cos_sigma = math.sin(sigma), math.cos(sigma)
    lat2 = math.atan2(sin_U1 * cos_sigma + cos_U1 * sin_sigma * math.cos(alpha1),
                      (1 - EARTH_F) * math.sqrt(
                          sin_alpha ** 2 + (sin_U1 * sin_sigma - cos_U1 * cos_sigma * math.cos(alpha1)) ** 2))
    lambda_rad = math.atan2(sin_sigma * math.sin(alpha1),
                            cos_U1 * cos_sigma - sin_U1 * sin_sigma * math.cos(alpha1))
    C = EARTH_F / 16 * cos2_alpha * (4 + EARTH_F * (4 - 3 * cos2_alpha))
    L = lambda_rad - (1 - C) * EARTH_F * sin_alpha * (sigma + C * sin_sigma *
                                                      (cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m ** 2)))
    return math.degrees(lat2), math.degrees(lon1 + L)


def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371000
    φ1, φ2 = math.radians(lat1), math.radians(lat2)
    Δφ, Δλ = math.radians(lat2 - lat1), math.radians(lon2 - lon1)
    a = math.sin(Δφ / 2) ** 2 + math.cos(φ1) * math.cos(φ2) * math.sin(Δλ / 2) ** 2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def rewind_trajectory(lat_end, lon_end, alt_end, heading, speed, vz, steps):
    lat, lon, alt = lat_end, lon_end, alt_end
    current_heading = heading
    dist_step = speed * KNOT_TO_MS
    alt_step = vz / 60.0

    for _ in range(steps):
        back_azimuth = (current_heading + 180) % 360
        prev_lat, prev_lon = direct_geodetic_problem(lat, lon, back_azimuth, dist_step)
        prev_alt = alt - alt_step
        new_heading = andoyer_azimuth(prev_lat, prev_lon, lat, lon)
        lat, lon, alt = prev_lat, prev_lon, prev_alt
        current_heading = new_heading
    return lat, lon, alt


def angle_difference(a1, a2):
    """Кратчайшая разница между углами в градусах"""
    diff = abs(a1 - a2) % 360
    return min(diff, 360 - diff)


def smooth_angle(current, target, factor):
    """Плавное изменение угла"""
    diff = angle_difference(current, target)
    if diff < 0.5:
        return target

    # Определяем направление поворота (кратчайший путь)
    delta = (target - current) % 360
    if delta > 180:
        delta = delta - 360

    # Плавный шаг
    step = delta * factor
    if abs(step) < 0.5:
        step = delta if abs(delta) < 0.5 else step

    result = current + step
    return result % 360


def get_safe_heading(current_heading, horizontal_bands):
    """Получает безопасное направление от DAIDALUS"""
    for band in horizontal_bands:
        if band.get('Bands_Type') == 'RECOVERY':
            center = (band['low'] + band['high']) / 2
            return center
    for band in horizontal_bands:
        if band.get('Bands_Type') == 'NONE':
            center = (band['low'] + band['high']) / 2
            return center
    return (current_heading + 90) % 360


def has_conflict_from_bands(horizontal_bands):
    if not horizontal_bands:
        return False
    for band in horizontal_bands:
        if band.get('Bands_Type') in ['NEAR', 'RECOVERY']:
            return True
    return False


def generate_single_daa_smooth(filename, ownship_init, ac1_init, duration_sec=60):
    set_scenario(filename)

    lat_o, lon_o, alt_o = ownship_init['lat'], ownship_init['lon'], ownship_init['alt']
    vx_o, vy_o = ownship_init['vx'], ownship_init['vy']
    vz_o_func = ownship_init['vz_func']

    lat_a, lon_a, alt_a = ac1_init['lat'], ac1_init['lon'], ac1_init['alt']
    vx_a, vy_a = ac1_init['vx'], ac1_init['vy']
    vz_a_func = ac1_init['vz_func']

    ownship_heading, ownship_speed = vx_vy_to_heading(vx_o, vy_o)
    intruder_heading, intruder_speed = vx_vy_to_heading(vx_a, vy_a)

    # СОХРАНЯЕМ ИСХОДНЫЙ КУРС (он не должен меняться!)
    ORIGINAL_HEADING = ownship_heading
    print(f"Исходный курс: {ORIGINAL_HEADING:.1f}°")

    # Состояние
    phase = "NORMAL"  # NORMAL, AVOIDANCE, RETURN
    target_heading = ORIGINAL_HEADING  # Целевой курс (изначально равен исходному)
    avoidance_heading = None  # Запомненный курс уклонения
    conflict_start_time = None
    return_start_time = None

    # Окно истории для определения линии пути
    history = deque(maxlen=20)
    original_positions = []  # Сохраняем позиции для возврата

    noise_std = 0.3
    NOISE_STD = {'lat_m': 5.0 * noise_std, 'lon_m': 5.0 * noise_std, 'alt_ft': 10.0 * noise_std,
                 'track_deg': 2.0 * noise_std, 'gs_knot': 2.0 * noise_std, 'vs_fpm': 5.0 * noise_std}

    with open(filename, 'w') as f:
        f.write("NAME, lat, lon, alt, vx, vy, vz, time\n")
        f.write("unitless, [deg], [deg], [ft], [knot], [knot], [fpm], [s]\n")

        for t in range(0, duration_sec + 1):
            vz_o = vz_o_func(t)
            vz_a = vz_a_func(t)

            f.write(f"Ownship, {lat_o:.8f}, {lon_o:.8f}, {alt_o:.0f}, "
                    f"{vx_o:.1f}, {vy_o:.1f}, {vz_o}, {t}\n")
            f.write(f"AC1, {lat_a:.8f}, {lon_a:.8f}, {alt_a:.0f}, "
                    f"{vx_a:.1f}, {vy_a:.1f}, {vz_a}, {t}\n")

            # Сохраняем позиции для истории
            history.append((lat_o, lon_o))
            original_positions.append((lat_o, lon_o, ownship_heading, t))

            # Расстояние до intruder
            horiz_dist = andoyer_distance(lat_o, lon_o, lat_a, lon_a)
            vert_dist = abs(alt_a - alt_o)

            # Получаем рекомендации от DAIDALUS
            lat_rad = math.radians(lat_o)
            std_lat_deg = NOISE_STD['lat_m'] / 111320.0
            std_lon_deg = NOISE_STD['lon_m'] / (111320.0 * math.cos(lat_rad))

            json_template = {
                "time": t,
                "ownship": {
                    "lat": lat_o + random.gauss(0, std_lat_deg),
                    "lon": lon_o + random.gauss(0, std_lon_deg),
                    "alt": alt_o + random.gauss(0, NOISE_STD['alt_ft']),
                    "track": (ownship_heading + random.gauss(0, NOISE_STD['track_deg'])) % 360,
                    "gs": max(0, ownship_speed + random.gauss(0, NOISE_STD['gs_knot'])),
                    "vs": vz_o + random.gauss(0, NOISE_STD['vs_fpm'])
                },
                "traffic": [{
                    "id": 'AC1',
                    "lat": lat_a + random.gauss(0, std_lat_deg),
                    "lon": lon_a + random.gauss(0, std_lon_deg),
                    "alt": alt_a + random.gauss(0, NOISE_STD['alt_ft']),
                    "track": (intruder_heading + random.gauss(0, NOISE_STD['track_deg'])) % 360,
                    "gs": max(0, intruder_speed + random.gauss(0, NOISE_STD['gs_knot'])),
                    "vs": vz_a + random.gauss(0, NOISE_STD['vs_fpm'])
                }]
            }

            horizontal_bands = run_simulation(json_template)
            daa_conflict = has_conflict_from_bands(horizontal_bands)
            real_conflict = horiz_dist < SAFE_HORIZ_M and vert_dist < SAFE_VERT_FT

            # ============================================================
            # ЛОГИКА СМЕНЫ ФАЗ
            # ============================================================

            if phase == "NORMAL":
                # Нормальный полёт - держим исходный курс
                target_heading = ORIGINAL_HEADING
                ownship_heading = smooth_angle(ownship_heading, ORIGINAL_HEADING, 0.1)

                if real_conflict or daa_conflict:
                    phase = "AVOIDANCE"
                    conflict_start_time = t
                    # Получаем безопасное направление
                    avoidance_heading = get_safe_heading(ownship_heading, horizontal_bands)
                    target_heading = avoidance_heading
                    print(f"\n=== t={t}: НАЧАЛО УКЛОНЕНИЯ ===")
                    print(f"  Расстояние: {horiz_dist:.0f}м")
                    print(f"  Исходный курс: {ORIGINAL_HEADING:.1f}° -> целевой: {target_heading:.1f}°")

            elif phase == "AVOIDANCE":
                # Уклонение - следуем безопасному курсу
                # Корректируем целевой курс при необходимости
                if horizontal_bands and (daa_conflict or horiz_dist < SAFE_HORIZ_M * 2):
                    new_target = get_safe_heading(ownship_heading, horizontal_bands)
                    if angle_difference(new_target, target_heading) > 15:
                        target_heading = new_target
                        print(f"t={t}: Обновлён целевой курс: {target_heading:.1f}°")

                # Плавно приближаем текущий курс к целевому
                ownship_heading = smooth_angle(ownship_heading, target_heading, HEADING_SMOOTHING)

                # Проверяем, можно ли вернуться (конфликт разрешён)
                if not daa_conflict and horiz_dist > SAFE_HORIZ_M * 3 and t - conflict_start_time > 5:
                    phase = "RETURN"
                    return_start_time = t
                    # Целевой курс для возврата - ИСХОДНЫЙ
                    target_heading = ORIGINAL_HEADING
                    print(f"\n=== t={t}: НАЧАЛО ВОЗВРАТА К ИСХОДНОМУ КУРСУ ===")
                    print(f"  Текущий курс: {ownship_heading:.1f}°, целевой: {target_heading:.1f}°")

            elif phase == "RETURN":
                # Возврат к исходному курсу
                target_heading = ORIGINAL_HEADING
                ownship_heading = smooth_angle(ownship_heading, ORIGINAL_HEADING, RETURN_SMOOTHING)

                # Завершаем возврат, когда курс достаточно близок к исходному
                if angle_difference(ownship_heading, ORIGINAL_HEADING) < 3:
                    phase = "NORMAL"
                    print(f"\n=== t={t}: ВОЗВРАТ ЗАВЕРШЁН ===")
                    print(f"  Итоговый курс: {ownship_heading:.1f}° (исходный: {ORIGINAL_HEADING:.1f}°)")

            # Вывод статуса
            status_map = {"NORMAL": "НОРМА", "AVOIDANCE": "УКЛОНЕНИЕ", "RETURN": "ВОЗВРАТ"}
            print(f"t={t:3d}: dist={horiz_dist:5.0f}м, {status_map[phase]:10}, курс={ownship_heading:5.1f}°")

            # ============================================================
            # ОБНОВЛЕНИЕ ПОЗИЦИЙ
            # ============================================================

            if t < duration_sec:
                distance_o = ownship_speed * KNOT_TO_MS
                new_lat_o, new_lon_o = direct_geodetic_problem(lat_o, lon_o, ownship_heading, distance_o)
                new_alt_o = alt_o + vz_o / 60.0
                new_ownship_heading = andoyer_azimuth(lat_o, lon_o, new_lat_o, new_lon_o)

                distance_a = intruder_speed * KNOT_TO_MS
                new_lat_a, new_lon_a = direct_geodetic_problem(lat_a, lon_a, intruder_heading, distance_a)
                new_alt_a = alt_a + vz_a / 60.0
                new_intruder_heading = andoyer_azimuth(lat_a, lon_a, new_lat_a, new_lon_a)

                lat_o, lon_o, alt_o = new_lat_o, new_lon_o, new_alt_o
                lat_a, lon_a, alt_a = new_lat_a, new_lon_a, new_alt_a
                ownship_heading = new_ownship_heading
                intruder_heading = new_intruder_heading

    # Проверка: совпадает ли начальный и конечный курс
    first_traj = original_positions[0] if original_positions else None
    last_traj = original_positions[-1] if original_positions else None

    if first_traj:
        print(f"\n=== ИТОГОВАЯ ПРОВЕРКА ===")
        print(f"Курс на 0 секунде: {first_traj[2]:.1f}°")
        print(f"Курс на {duration_sec} секунде: {ownship_heading:.1f}°")

        if angle_difference(first_traj[2], ownship_heading) < 5:
            print("КУРС СОХРАНЁН!")
        else:
            print(f"ВНИМАНИЕ: Курс изменился на {angle_difference(first_traj[2], ownship_heading):.1f}°")

    finalize_scenario()
    print(f"\nГенерация завершена")


def generate_conflict_scenario_smooth(filename, duration_sec=60, conflict_time=20,
                                      region=None, ownship_config=None, ac1_config=None):
    if region is None:
        region = {'lat_center': 45.0, 'lat_spread': 0.01, 'lon_center': 45.0,
                  'lon_spread': 0.02, 'alt_center': 5000, 'alt_spread': 1000}

    if ownship_config is None:
        ownship_config = Ownship(speed_knot=120, altitude=5000, heading_deg=90, vz=0)
    if ac1_config is None:
        ac1_config = Intruder(speed_knot=150, heading_deg=270, vz=0)

    # Позиция в момент конфликта
    lat_c = random.uniform(region['lat_center'] - region['lat_spread'],
                           region['lat_center'] + region['lat_spread'])
    lon_c = random.uniform(region['lon_center'] - region['lon_spread'],
                           region['lon_center'] + region['lon_spread'])
    alt_c = ownship_config.altitude

    # Позиция AC1
    r = random.uniform(0, SAFE_HORIZ_M)
    angle = random.uniform(0, 2 * math.pi)
    dx_m, dy_m = r * math.cos(angle), r * math.sin(angle)

    lat_c_rad = math.radians(lat_c)
    sin_lat = math.sin(lat_c_rad)
    ro1 = (EARTH_A * (1 - EARTH_E2)) / (math.sqrt((1 - EARTH_E2 * sin_lat ** 2) ** 3)) + alt_c
    ro2 = EARTH_A / math.sqrt(1 - EARTH_E2 * sin_lat ** 2) + alt_c

    dlat = math.degrees(dy_m / ro1)
    ac_lat_c = lat_c + dlat
    dlon = dx_m / (ro2 * math.cos(math.radians(ac_lat_c)))
    ac_lon_c = lon_c + math.degrees(dlon)
    ac_alt_c = alt_c

    # Параметры движения
    own_heading = ownship_config.heading_deg
    own_speed = ownship_config.speed_knot
    own_vz = ownship_config.vz
    own_vx, own_vy = heading_to_vx_vy(own_heading, own_speed)

    ac_heading = ac1_config.heading_deg
    ac_speed = ac1_config.speed_knot
    ac_vz = ac1_config.vz
    ac_vx, ac_vy = heading_to_vx_vy(ac_heading, ac_speed)

    # Отмотка назад
    lat_o0, lon_o0, alt_o0 = rewind_trajectory(lat_c, lon_c, alt_c, own_heading, own_speed, own_vz, conflict_time)
    lat_a0, lon_a0, alt_a0 = rewind_trajectory(ac_lat_c, ac_lon_c, ac_alt_c, ac_heading, ac_speed, ac_vz, conflict_time)

    def make_vz_func(vz_val):
        return lambda t: int(round(vz_val))

    ownship_init = {'lat': lat_o0, 'lon': lon_o0, 'alt': alt_o0,
                    'vx': own_vx, 'vy': own_vy, 'vz_func': make_vz_func(own_vz)}
    ac1_init = {'lat': lat_a0, 'lon': lon_a0, 'alt': alt_a0,
                'vx': ac_vx, 'vy': ac_vy, 'vz_func': make_vz_func(ac_vz)}

    generate_single_daa_smooth(filename, ownship_init, ac1_init, duration_sec)


def generate_multiple_scenarios(output_dir="smooth_conflict_scenarios", num_scenarios=5,
                                duration_sec=60, conflict_time=20, seed=2025,
                                region=None, ownship_config=None, ac1_config=None):
    if seed is not None:
        random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    for i in range(1, num_scenarios + 1):
        fname = os.path.join(output_dir, f"conflict_{i:04d}.daa")
        generate_conflict_scenario_smooth(fname, duration_sec, conflict_time,
                                          region, ownship_config, ac1_config)
        print(f"Создан: {fname}")


if __name__ == "__main__":
    ownship_params = Ownship(heading_deg=90, speed_knot=120, altitude=5000, vz=0)
    ac1_params = Intruder(heading_deg=270, speed_knot=150, vz=0)

    region = {'lat_center': 45.0, 'lat_spread': 0.01, 'lon_center': 45.0,
              'lon_spread': 0.02, 'alt_center': 5000, 'alt_spread': 1000}

    generate_multiple_scenarios(
        output_dir="smooth_conflict_scenarios",
        num_scenarios=3,
        duration_sec=100,
        conflict_time=20,
        seed=2025,
        region=region,
        ownship_config=ownship_params,
        ac1_config=ac1_params
    )