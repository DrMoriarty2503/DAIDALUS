import requests
import time
import json

URL = "http://localhost:8080/detect"
HEALTH_URL = "http://localhost:8080/health"


def run_simulation(payload):
    current_time = time.time()
    payload_json = json.dumps(payload)
    try:
        r = requests.get(HEALTH_URL, timeout=2)
        print(f"Server: {r.text}")
    except Exception as e:
        print(f"Connection failed: {e}")
        return

    try:
            resp = requests.post(URL, json=payload, timeout=1)

            if resp.status_code == 200:
                data = resp.json()

                print(f"STAP {payload['time']}")

                # Alerting Logic
                alert = data.get("Alerting Logic", "Unknown")
                print(f"Alerting Logic: {alert}")

                print_bands(data.get("Horizontal Bands", []), "Horizontal Bands")
                print_bands(data.get("Horizontal Speed Bands", []), "Horizontal Speed Bands")
                print_bands(data.get("Vertical Bands", []), "Vertical Bands")
                print_bands(data.get("Vertical Speed Bands", []), "Vertical Speed Bands")

                print()

            else:
                print(f"Tick {payload['time']}: Error {resp.status_code}")


            time.sleep(1)

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