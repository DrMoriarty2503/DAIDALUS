import requests
import time
import json
import csv
import os

URL = "http://localhost:8080/detect"
HEALTH_URL = "http://localhost:8080/health"

# Глобальная переменная для хранения данных текущего сценария
current_scenario_buffer = []
current_scenario_name = None


def set_scenario(scenario_name):
    """Устанавливает новый сценарий и очищает буфер"""
    global current_scenario_buffer, current_scenario_name
    current_scenario_buffer = []
    current_scenario_name = scenario_name


def save_scenario_to_csv():
    """Сохраняет собранные данные текущего сценария в CSV"""
    global current_scenario_buffer, current_scenario_name

    if not current_scenario_buffer:
        print(f"Нет данных для сохранения: {current_scenario_name}")
        return

    if current_scenario_name is None:
        filename = "recommendations_unknown.csv"
    else:
        # Заменяем .daa на .csv
        filename = current_scenario_name.replace('.daa', '_recommendations.csv')

    # Определяем все возможные колонки
    all_keys = set()
    for row in current_scenario_buffer:
        all_keys.update(row.keys())

    sorted_keys = sorted(all_keys)

    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=sorted_keys)
        writer.writeheader()
        writer.writerows(current_scenario_buffer)

    print(f"✅ Сохранены рекомендации: {filename} ({len(current_scenario_buffer)} записей)")


def run_simulation(payload):
    """Отправляет данные в DAIDALUS и собирает рекомендации"""
    global current_scenario_buffer

    try:
        # Проверка здоровья сервера (только один раз)
        if len(current_scenario_buffer) == 0:
            try:
                r = requests.get(HEALTH_URL, timeout=2)
                print(f"Server: {r.text}")
            except Exception as e:
                print(f"Connection failed: {e}")
                return

        resp = requests.post(URL, json=payload, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            # print(data)
            # print(data)
            current_time = payload.get('time', 0)

            print(f"STAP {current_time}")

            # Alerting Logic
            alert = data.get("Alerting Logic", "Unknown")
            print(f"Alerting Logic: {alert}")

            # Извлекаем bands
            horizontal_bands = data.get("Horizontal Bands", [])
            horizontal_speed_bands = data.get("Horizontal Speed Bands", [])
            vertical_bands = data.get("Vertical Bands", [])
            vertical_speed_bands = data.get("Vertical Speed Bands", [])

            print_bands(horizontal_bands, "Horizontal Bands")
            print_bands(horizontal_speed_bands, "Horizontal Speed Bands")
            print_bands(vertical_bands, "Vertical Bands")
            print_bands(vertical_speed_bands, "Vertical Speed Bands")
            # Формируем запись для CSV
            csv_row = {
                'time': current_time,
                'alert_level': alert,
                'ownship_lat': payload['ownship']['lat'],
                'ownship_lon': payload['ownship']['lon'],
                'ownship_alt': payload['ownship']['alt'],
                'ownship_track': payload['ownship']['track'],
                'ownship_gs': payload['ownship']['gs'],
                'ownship_vs': payload['ownship']['vs'],
            }

            # Добавляем информацию о трафике
            if payload.get('traffic') and len(payload['traffic']) > 0:
                traffic = payload['traffic'][0]
                csv_row['traffic_id'] = traffic.get('id', '')
                csv_row['traffic_lat'] = traffic.get('lat')
                csv_row['traffic_lon'] = traffic.get('lon')
                csv_row['traffic_alt'] = traffic.get('alt')
                csv_row['traffic_track'] = traffic.get('track')
                csv_row['traffic_gs'] = traffic.get('gs')
                csv_row['traffic_vs'] = traffic.get('vs')

            # Добавляем Horizontal Bands
            for i, band in enumerate(horizontal_bands):
                band_type = band.get("Bands_Type", "UNKNOWN")
                csv_row[f'h_band_{i}_type'] = band_type
                csv_row[f'h_band_{i}_low'] = band.get("low")
                csv_row[f'h_band_{i}_high'] = band.get("high")
                csv_row[f'h_band_{i}_unit'] = band.get("unit", "")

            # Добавляем Horizontal Speed Bands
            for i, band in enumerate(horizontal_speed_bands):
                band_type = band.get("Bands_Type", "UNKNOWN")
                csv_row[f'hs_band_{i}_type'] = band_type
                csv_row[f'hs_band_{i}_low'] = band.get("low")
                csv_row[f'hs_band_{i}_high'] = band.get("high")
                csv_row[f'hs_band_{i}_unit'] = band.get("unit", "")

            # Добавляем Vertical Bands
            for i, band in enumerate(vertical_bands):
                band_type = band.get("Bands_Type", "UNKNOWN")
                csv_row[f'v_band_{i}_type'] = band_type
                csv_row[f'v_band_{i}_low'] = band.get("low")
                csv_row[f'v_band_{i}_high'] = band.get("high")
                csv_row[f'v_band_{i}_unit'] = band.get("unit", "")

            # Добавляем Vertical Speed Bands
            for i, band in enumerate(vertical_speed_bands):
                band_type = band.get("Bands_Type", "UNKNOWN")
                csv_row[f'vs_band_{i}_type'] = band_type
                csv_row[f'vs_band_{i}_low'] = band.get("low")
                csv_row[f'vs_band_{i}_high'] = band.get("high")
                csv_row[f'vs_band_{i}_unit'] = band.get("unit", "")

            current_scenario_buffer.append(csv_row)

        else:
            print(f"Tick {payload.get('time', '?')}: Error {resp.status_code}")

    except requests.exceptions.Timeout:
        print(f"Timeout на t={payload.get('time', '?')}, пропускаем")
    except requests.exceptions.ConnectionError:
        print("Сервер недоступен, проверьте Docker контейнер")
    except KeyboardInterrupt:
        print("\nStopped")
    except Exception as e:
        print(f"Request error: {e}")


    return horizontal_bands


def finalize_scenario():
    """Сохраняет данные текущего сценария и очищает буфер"""
    save_scenario_to_csv()
    global current_scenario_buffer
    current_scenario_buffer = []


def print_bands(bands, title):
    """Выводит bands в консоль"""
    print(f"{title}:")
    if not bands:
        print("  (no data)")
        return
    for band in bands:
        b_type = band.get("Bands_Type", "UNKNOWN")
        low = band.get("low", "?")
        high = band.get("high", "?")
        unit = band.get("unit", "")
        if unit:
            print(f"  {b_type}: {low}-{high} {unit}")
        else:
            print(f"  {b_type}: {low}-{high}")


if __name__ == "__main__":
    run_simulation()