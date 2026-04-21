import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math


def read_daa_trajectory(filename):
    """Читает .daa файл и возвращает траектории"""
    ownship_lats = []
    ownship_lons = []
    intruder_lats = []
    intruder_lons = []
    times = []

    with open(filename, 'r', encoding='utf-8') as f:
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
            t = float(parts[7])
        except ValueError:
            continue

        if name == 'Ownship':
            ownship_lats.append(lat)
            ownship_lons.append(lon)
            if len(times) <= len(ownship_lats):
                times.append(t)
        elif name == 'AC1':
            intruder_lats.append(lat)
            intruder_lons.append(lon)

    return {
        'ownship_lats': ownship_lats,
        'ownship_lons': ownship_lons,
        'intruder_lats': intruder_lats,
        'intruder_lons': intruder_lons,
        'times': times
    }


def calculate_distance(lat1, lon1, lat2, lon2):
    """Расстояние между точками в метрах"""
    R = 6371000
    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)

    a = math.sin(dlat / 2) ** 2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


def animate_trajectory(filename, speed=100.0, save_gif=False):
    """
    Анимирует движение двух ВС.

    Параметры:
        filename - путь к .daa файлу
        speed - скорость анимации (1.0 = реальное время, 2.0 = в 2 раза быстрее)
        save_gif - сохранить ли GIF файл
    """

    # Загружаем данные
    data = read_daa_trajectory(filename)
    ownship_lats = data['ownship_lats']
    ownship_lons = data['ownship_lons']
    intruder_lats = data['intruder_lats']
    intruder_lons = data['intruder_lons']
    times = data['times']

    # Создаем фигуру с двумя графиками
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # ============================================================
    # ЛЕВЫЙ ГРАФИК: Анимация движения (вид сверху)
    # ============================================================

    # Находим границы для масштабирования
    all_lats = ownship_lats + intruder_lats
    all_lons = ownship_lons + intruder_lons
    lat_center = (max(all_lats) + min(all_lats)) / 2
    lon_center = (max(all_lons) + min(all_lons)) / 2
    lat_range = max(all_lats) - min(all_lats)
    lon_range = max(all_lons) - min(all_lons)
    max_range = max(lat_range, lon_range) * 0.6

    ax1.set_xlim(lon_center - max_range, lon_center + max_range)
    ax1.set_ylim(lat_center - max_range, lat_center + max_range)

    # Траектории (следы)
    trail_ownship, = ax1.plot([], [], 'b-', linewidth=1, alpha=0.5, label='Траектория Ownship')
    trail_intruder, = ax1.plot([], [], 'r-', linewidth=1, alpha=0.5, label='Траектория Intruder')

    # Точки (самолеты)
    ownship_point, = ax1.plot([], [], 'bo', markersize=10, label='Ownship')
    intruder_point, = ax1.plot([], [], 'rs', markersize=10, label='Intruder (AC1)')

    # Линия расстояния между ВС
    distance_line, = ax1.plot([], [], 'g--', linewidth=1, alpha=0.7)

    # Текст с информацией
    info_text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, fontsize=10,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax1.set_xlabel('Долгота (градусы)')
    ax1.set_ylabel('Широта (градусы)')
    ax1.set_title('Движение воздушных судов (вид сверху)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')

    # ============================================================
    # ПРАВЫЙ ГРАФИК: Расстояние по времени
    # ============================================================

    # Рассчитываем все расстояния
    distances = []
    for i in range(len(times)):
        dist = calculate_distance(ownship_lats[i], ownship_lons[i],
                                  intruder_lats[i], intruder_lons[i])
        distances.append(dist)

    ax2.plot(times, distances, 'g-', linewidth=2, alpha=0.7)
    ax2.axhline(y=150, color='r', linestyle='--', linewidth=1)
    ax2.axhline(y=609.9, color='orange', linestyle='--', linewidth=1)

    # Текущая точка на графике расстояния
    current_point, = ax2.plot([], [], 'ro', markersize=8)
    distance_text = ax2.text(0.02, 0.95, '', transform=ax2.transAxes, fontsize=10,
                             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax2.set_xlabel('Время (секунды)')
    ax2.set_ylabel('Расстояние (метры)')
    ax2.set_title('Расстояние между ВС')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, max(times))
    ax2.set_ylim(0, max(distances) * 1.1)

    # ============================================================
    # Функция анимации
    # ============================================================

    # Интервал между кадрами в миллисекундах (1 секунда реального времени / speed)
    interval = 1000 / speed

    def update(frame):
        """Обновляет кадр анимации"""

        # Текущие координаты (берем до текущего кадра)
        current_own_lats = ownship_lats[:frame + 1]
        current_own_lons = ownship_lons[:frame + 1]
        current_int_lats = intruder_lats[:frame + 1]
        current_int_lons = intruder_lons[:frame + 1]

        # Обновляем следы (траектории)
        trail_ownship.set_data(current_own_lons, current_own_lats)
        trail_intruder.set_data(current_int_lons, current_int_lats)

        # Обновляем позиции самолетов
        ownship_point.set_data([ownship_lons[frame]], [ownship_lats[frame]])
        intruder_point.set_data([intruder_lons[frame]], [intruder_lats[frame]])

        # Обновляем линию расстояния
        distance_line.set_data([ownship_lons[frame], intruder_lons[frame]],
                               [ownship_lats[frame], intruder_lats[frame]])

        # Рассчитываем текущее расстояние
        current_dist = calculate_distance(ownship_lats[frame], ownship_lons[frame],
                                          intruder_lats[frame], intruder_lons[frame])

        # Обновляем текст с информацией
        info_text.set_text(
            f'Время: {times[frame]:.0f} с\n'
            f'Ownship: ({ownship_lats[frame]:.4f}°, {ownship_lons[frame]:.4f}°)\n'
            f'Intruder: ({intruder_lats[frame]:.4f}°, {intruder_lons[frame]:.4f}°)\n'
            f'Расстояние: {current_dist:.1f} м'
        )

        # Обновляем график расстояния
        current_point.set_data([times[frame]], [distances[frame]])

        # Цвет текста в зависимости от опасности
        if current_dist < 150:
            color = 'red'
            status = '⚠️ КОНФЛИКТ NMAC!'
        elif current_dist < 609.9:
            color = 'orange'
            status = '⚠️ LOWC!'
        else:
            color = 'green'
            status = '✅ Безопасно'

        distance_text.set_text(f'Текущее расстояние: {current_dist:.1f} м\n{status}')
        distance_text.set_color(color)

        # Меняем цвет точек при конфликте
        if current_dist < 150:
            ownship_point.set_color('red')
            intruder_point.set_color('red')
        elif current_dist < 609.9:
            ownship_point.set_color('orange')
            intruder_point.set_color('orange')
        else:
            ownship_point.set_color('blue')
            intruder_point.set_color('red')

        return trail_ownship, trail_intruder, ownship_point, intruder_point, distance_line, info_text, current_point, distance_text

    # Создаем анимацию
    anim = animation.FuncAnimation(fig, update, frames=len(times),
                                   interval=interval, repeat=True, blit=False)

    plt.tight_layout()

    # Сохраняем GIF если нужно
    if save_gif:
        output_file = filename.replace('.daa', '_animation.gif')
        anim.save(output_file, writer='pillow', fps=10)
        print(f"Анимация сохранена: {output_file}")

    plt.show()
    return anim


def interactive_timeline(filename):
    """
    Интерактивная версия с ползунком времени.
    Позволяет вручную перемещаться по времени.
    """
    import matplotlib.widgets as widgets

    data = read_daa_trajectory(filename)
    ownship_lats = data['ownship_lats']
    ownship_lons = data['ownship_lons']
    intruder_lats = data['intruder_lats']
    intruder_lons = data['intruder_lons']
    times = data['times']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Настройка левого графика
    all_lats = ownship_lats + intruder_lats
    all_lons = ownship_lons + intruder_lons
    lat_center = (max(all_lats) + min(all_lats)) / 2
    lon_center = (max(all_lons) + min(all_lons)) / 2
    lat_range = max(all_lats) - min(all_lats)
    lon_range = max(all_lons) - min(all_lons)
    max_range = max(lat_range, lon_range) * 0.6

    ax1.set_xlim(lon_center - max_range, lon_center + max_range)
    ax1.set_ylim(lat_center - max_range, lat_center + max_range)

    # Траектории
    ax1.plot(ownship_lons, ownship_lats, 'b-', linewidth=1, alpha=0.3)
    ax1.plot(intruder_lons, intruder_lats, 'r-', linewidth=1, alpha=0.3)

    # Точки (изначально на t=0)
    ownship_point, = ax1.plot(ownship_lons[0], ownship_lats[0], 'bo', markersize=10)
    intruder_point, = ax1.plot(intruder_lons[0], intruder_lats[0], 'rs', markersize=10)
    distance_line, = ax1.plot([ownship_lons[0], intruder_lons[0]],
                              [ownship_lats[0], intruder_lats[0]], 'g--', linewidth=1)

    info_text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, fontsize=10,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax1.set_xlabel('Долгота (градусы)')
    ax1.set_ylabel('Широта (градусы)')
    ax1.set_title('Движение ВС (перемещайте ползунок)')
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')

    # Правый график
    distances = [calculate_distance(ownship_lats[i], ownship_lons[i],
                                    intruder_lats[i], intruder_lons[i])
                 for i in range(len(times))]

    ax2.plot(times, distances, 'g-', linewidth=2)
    ax2.axhline(y=150, color='r', linestyle='--', label='NMAC (150 м)')
    ax2.axhline(y=609.9, color='orange', linestyle='--', label='LOWC (610 м)')
    ax2.set_xlabel('Время (секунды)')
    ax2.set_ylabel('Расстояние (метры)')
    ax2.set_title('Расстояние между ВС')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Текущая точка на графике расстояния
    current_point, = ax2.plot(times[0], distances[0], 'ro', markersize=8)
    distance_text = ax2.text(0.02, 0.95, '', transform=ax2.transAxes, fontsize=10,
                             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Ползунок времени
    ax_slider = plt.axes([0.2, 0.02, 0.6, 0.03])
    time_slider = widgets.Slider(ax_slider, 'Время (с)', 0, len(times) - 1, valinit=0, valstep=1)

    def update(val):
        frame = int(time_slider.val)

        ownship_point.set_data([ownship_lons[frame]], [ownship_lats[frame]])
        intruder_point.set_data([intruder_lons[frame]], [intruder_lats[frame]])
        distance_line.set_data([ownship_lons[frame], intruder_lons[frame]],
                               [ownship_lats[frame], intruder_lats[frame]])

        current_dist = distances[frame]

        info_text.set_text(
            f'Время: {times[frame]:.0f} с\n'
            f'Ownship: ({ownship_lats[frame]:.4f}°, {ownship_lons[frame]:.4f}°)\n'
            f'Intruder: ({intruder_lats[frame]:.4f}°, {intruder_lons[frame]:.4f}°)\n'
            f'Расстояние: {current_dist:.1f} м'
        )

        current_point.set_data([times[frame]], [current_dist])

        if current_dist < 150:
            color = 'red'
            status = '⚠️ КОНФЛИКТ NMAC!'
        elif current_dist < 609.9:
            color = 'orange'
            status = '⚠️ LOWC!'
        else:
            color = 'green'
            status = '✅ Безопасно'

        distance_text.set_text(f'Расстояние: {current_dist:.1f} м\n{status}')
        distance_text.set_color(color)

        fig.canvas.draw_idle()

    time_slider.on_changed(update)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Использование:")
        print("  python animate.py <файл.daa> [--speed N] [--gif] [--interactive]")
        print("\nПримеры:")
        print("  python animate.py random_conflicts/conflict_0001.daa")
        print("  python animate.py random_conflicts/conflict_0001.daa --speed 2.0")
        print("  python animate.py random_conflicts/conflict_0001.daa --gif")
        print("  python animate.py random_conflicts/conflict_0001.daa --interactive")
        sys.exit(1)

    filename = sys.argv[1]
    speed = 1.0
    save_gif = '--gif' in sys.argv
    interactive = '--interactive' in sys.argv

    # Поиск параметра скорости
    for i, arg in enumerate(sys.argv):
        if arg == '--speed' and i + 1 < len(sys.argv):
            try:
                speed = float(sys.argv[i + 1])
            except ValueError:
                pass

    if interactive:
        interactive_timeline(filename)
    else:
        animate_trajectory(filename, speed=speed, save_gif=save_gif)