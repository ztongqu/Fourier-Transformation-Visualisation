import FAE4D333_maths, math
import matplotlib.pyplot as plt
import threading, time

start_time = time.time()
func: type = type(print)

def Lerp(A: complex, B: complex, t: float) -> complex:
    return (1 - t) * A + t * B
def Lerp2(A: complex, B: complex, C: complex, t: float) -> complex:
    return Lerp(Lerp(A, B, t), Lerp(B, C, t), t)
def Lerp3(A: complex, B: complex, C: complex, D: complex, t: float) -> complex:
    return Lerp(Lerp2(A, B, C, t), Lerp2(B, C, D, t), t)

def line(x: float) -> complex:
    return x
def I(t: float) -> complex:
    ps: list[complex] = [-3+8j, 3+8j, 3+6j, 1+6j, 1-6j, 3-6j, 3-8j, -3-8j, -3-6j, -1-6j, -1+6j, -3+6j]
    if t < 1/12:
        return Lerp(ps[0], ps[1], 12*t)
    elif t < 2/12:
        return Lerp(ps[1], ps[2], 12*t-1)
    elif t < 3/12:
        return Lerp(ps[2], ps[3], 12*t-2)
    elif t < 4/12:
        return Lerp(ps[3], ps[4], 12*t-3)
    elif t < 5/12:
        return Lerp(ps[4], ps[5], 12*t-4)
    elif t < 6/12:
        return Lerp(ps[5], ps[6], 12*t-5)
    elif t < 7/12:
        return Lerp(ps[6], ps[7], 12*t-6)
    elif t < 8/12:
        return Lerp(ps[7], ps[8], 12*t-7)
    elif t < 9/12:
        return Lerp(ps[8], ps[9], 12*t-8)
    elif t < 10/12:
        return Lerp(ps[9], ps[10], 12*t-9)
    elif t < 11/12:
        return Lerp(ps[10], ps[11], 12*t-10)
    else:
        return Lerp(ps[11], ps[0], 12*t-11)
def R(t: float) -> complex:
    ps: list[complex] = [-4 + 8j, 8j, 8 + 8j, 8 + 1j, 2, 6 - 8j, 4 - 8j, 0, -2, -2 - 8j, -4 - 8j, -2 + 6j, -2 + 2j, 1 + 2j, 5 + 2j, 5 + 6j, 1 + 6j]
    if t < 1 / 13:
        return Lerp(ps[0], ps[1], 13 * t)
    elif t < 2 / 13:
        return Lerp3(ps[1], ps[2], ps[3], ps[4], 13 * t - 1)
    elif t < 3 / 13:
        return Lerp(ps[4], ps[5], 13 * t - 2)
    elif t < 4 / 13:
        return Lerp(ps[5], ps[6], 13 * t - 3)
    elif t < 5 / 13:
        return Lerp(ps[6], ps[7], 13 * t - 4)
    elif t < 6 / 13:
        return Lerp(ps[7], ps[8], 13 * t - 5)
    elif t < 7 / 13:
        return Lerp(ps[8], ps[9], 13 * t - 6)
    elif t < 8 / 13:
        return Lerp(ps[9], ps[10], 13 * t - 7)
    elif t < 9 / 13:
        return Lerp(ps[10], ps[0], 13 * t - 8)
    elif t < 10 / 13:
        return Lerp(ps[11], ps[12], 13 * t - 9)
    elif t < 11 / 13:
        return Lerp(ps[12], ps[13], 13 * t - 10)
    elif t < 12 / 13:
        return Lerp3(ps[13], ps[14], ps[15], ps[16], 13 * t - 11)
    else:
        return Lerp(ps[16], ps[11], 13 * t - 12)

skip = 1
if not bool(skip):
    x_points: list[float] = []
    y_points: list[float] = []
    for i in range(0, 1001):
        x_points.append(R(i/1000).real)
        y_points.append(R(i/1000).imag)
    plt.plot(x_points, y_points)
    plt.grid()
    plt.axis("equal")
    plt.show()


class FT_Thread(threading.Thread):
    def __init__(self, thread_id: int):
        threading.Thread.__init__(self)
        self.thread_id = thread_id

    def run(self):
        global coef, dic
        step: int = int(dic["range"] * 2 / dic["thread_num"])
        c_value: int = (2 * dic["range"]) % step
        start: int = step * self.thread_id - dic["range"]
        end: int = min(step + start, dic["range"])
        print(f"thread {self.thread_id} is integrating from {start} to {end - 1}")
        for _ in range(start, end):
            if _ & 1 == 9:
                r = 0
            else:
                r = FAE4D333_maths.fourier_transform(dic["f"], _, dic["dx"], dic["acc"])
            coef[_] = r
"""
def riemann_sum(f: func, start: float, end: float, dx: float, acc: int | None = None, power: int | None = 0, thread_id: int | None = None) -> complex:
    global thread_num
    if thread_id is None:
        thread_id = 0
        thread_num = 1
    r = end - start
    _ = start + r / thread_num * thread_id
    s = 0
    # print("integrating from", _, "to", start + r / thread_num * (thread_id + 1))
    while _ < start + r / thread_num * (thread_id + 1):
        s += 0.5 * dx * (f(_) * cmath.exp(2j * cmath.pi * power * _) + f(_ + dx) * cmath.exp(2j * cmath.pi * power * (_ + dx)))
        _ += dx
    if acc is None:
        return s
    else:
        return complex(round(s.real, acc), round(s.imag, acc))
"""
dic = {"f": I, "dx": 1 / 200000, "acc": 4, "range": 50, "thread_num": 8}
coef: dict[int, complex] = {}

integration_threads = []
for i in range(0, dic["thread_num"]):
    integration_threads.append(FT_Thread(i))
    integration_threads[i].start()
for i in range(0, dic["thread_num"]):
    integration_threads[i].join()

end_time = time.time()
print("time:", end_time - start_time)
print(coef)
"""
start_time = time.time()
for i in range(-50, 50):
    if i & 1 == 9:
        coef[i] = 0
    else:
        coef[i] = FAE4D333_maths.fourier_transform(R, -i, 1/1000000, 5)
end_time = time.time()
print(end_time - start_time)
print(coef)
"""


