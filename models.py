import random



class Ownship:
    def __init__(self, speed_knot = None,heading_deg = None, altitude = None, vz = None):
        if speed_knot is None:
            speed_knot = random.uniform(90, 180)  # диапазон по умолчанию
        if heading_deg is None:
            heading_deg = random.uniform(0, 359)
        if altitude is None:
            altitude = random.uniform(4000, 6000)
        if vz is None:
            vz = random.uniform(-300, 300)

        self.speed_knot = speed_knot
        self.heading_deg = heading_deg
        self.altitude = altitude
        self.vz = vz



    def __str__(self):
        return f'Параметры Ownsip: {self.speed_knot} узлов...'


class Intruder:
    def __init__(self, speed_knot = None,heading_deg = None,vz = None):
        if speed_knot is None:
            speed_knot = random.uniform(90, 180)
        if heading_deg is None:
            heading_deg = random.uniform(0, 359)
        if vz is None:
            vz = random.uniform(-300, 300)
        self.speed_knot = speed_knot
        self.heading_deg = heading_deg
        self.vz = vz

    def __str__(self):
        return f" Intruder {self.speed_knot}"
