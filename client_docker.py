import requests
import time
import json

URL = "http://localhost:8080/detect"
HEALTH_URL = "http://localhost:8080/health"


def run_simulation():
    current_time = time.time()

    tick_rate = 0.5  # 2 Гц
    MAX_DURATION = 20  # Длительность в секундах
    MAX_TICKS = int(MAX_DURATION / tick_rate)  # = 40 тактов

    # Начальные позиции
    ownship = {
        "lat": 44.0, "lon": 44.0, "alt": 100.0,
        "track": 60.0, "gs": 50.0, "vs": 0.0
    }
    traffic = [{
        "id": "TRF1",
        "lat": 44.01, "lon": 45.01, "alt": 100.0,
        "track": 306.0, "gs": 50.0, "vs": 0.0
    }]

    # Проверка связи
    try:
        r = requests.get(HEALTH_URL, timeout=2)
        print(f"Server: {r.text}")
    except Exception as e:
        print(f"Connection failed: {e}")
        return

    try:

        tick = 0
        for tick in range(MAX_TICKS):
            payload = {
                "time": 0,
                "ownship": {
                    "lat": 45.0, "lon": 45.0, "alt": 5000,
                    "track": 0, "gs": 150, "vs": 0
                },
                "traffic": [{
                    "id": "AC1",
                    "lat": 45.0045, "lon": 45.0, "alt": 5000,
                    "track": 180, "gs": 150, "vs": 0
                }]
            }

            resp = requests.post(URL, json=payload, timeout=1)

            if resp.status_code == 200:
                data = resp.json()

                print(f"STAP {tick}")

                # Alerting Logic
                alert = data.get("Alerting Logic", "Unknown")
                print(f"Alerting Logic: {alert}")

                print_bands(data.get("Horizontal Bands", []), "Horizontal Bands")
                print_bands(data.get("Horizontal Speed Bands", []), "Horizontal Speed Bands")
                print_bands(data.get("Vertical Bands", []), "Vertical Bands")
                print_bands(data.get("Vertical Speed Bands", []), "Vertical Speed Bands")

                print()

            else:
                print(f"Tick {tick}: Error {resp.status_code}")

            ownship["lon"] += 0.0005
            traffic[0]["lon"] -= 0.0005
            current_time += tick_rate
            tick += 1

            time.sleep(tick_rate)

    except KeyboardInterrupt:
        print("\nStopped")


def print_bands(bands, title):
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