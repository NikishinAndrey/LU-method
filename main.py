import random as r
import numpy as np
import matplotlib.pyplot as plt

length = 10  # Размерность матрицы
length3 = 6  # Количество коэффициентов отделимости чисел
length4 = 10  # Количество дефформируемых матриц



# Создаем диагональную матрицу


def do_matrix(koef):
    mat_c = [[0 for i in range(length)] for j in range(length)]
    for i in range(length):
        mat_c[i][i] = (i * koef + 1)
        # mat_c[i][i] = 1 - i * koef
    mat_c = np.array(mat_c)  # Диагональная матрица
    a__ = np.array([[r.randint(1, 100) for i in range(length)] for j in range(length)])  # Создаем рандомную матрицу
    a_ = np.array(np.linalg.inv(a__))  # Находим обратную матрицу
    matrix = a__.dot(mat_c).dot(a_)  # Получаем конечную симметричную матрицу
    return matrix

# Записываем результат в переменную
# Диагональная матрица
# Собственные числа матрицы
# Матрица с точным значением собственных векторов
# print('Заданые собственные числа: ')
# print(eig_value1)
# print('Заданные собственные вектора: ')
# print(eig_vector1)
# a__ = np.array([[r.randint(1, 100) for i in range(length)] for j in range(length)])  # Создаем рандомную матрицу
# a_ = np.array(np.linalg.inv(a__))  # Находим обратную матрицу
# matrix = a__.dot(mat_c).dot(a_)  # Получаем конечную симметричную матрицу


def deforming_matrix(a):
    massive_a = [[[0 for i in range(length)] for j in range(length)] for k in range(length4)]
    for k in range(length4):
        for i in range(length):
            for j in range(length):
                massive_a[k][i][j] = a[i][j] * (1 + 0.01 * k)
    return massive_a

# Реализуем lu-разложение(хотя можно было воспользоваться встроенной функцией)


def lu_decomposition(a):
    u = [[0 for i in range(length)] for j in range(length)]
    l = [[0 for i in range(length)] for j in range(length)]
    for i in range(length):
        for j in range(i, length):
            sum1 = 0
            for k in range(i):
                sum1 += l[i][k] * u[k][j]
            u[i][j] = a[i][j] - sum1
            sum2 = 0
            for k in range(j):
                sum2 += l[j][k] * u[k][i]
            l[j][i] = (a[j][i] - sum2) / u[i][i]
    return l, u


# Реализуем lu-алгоритм нахождения собственных чисел с заданной точностью и максимальным количеством итераций


def lu_algorithm(a, e, n_max):
    current = [0 for j in range(length)]
    previous = [0 for k in range(length)]
    iteration = 0
    cur = 1
    while cur > e and iteration < n_max:
        res = lu_decomposition(a)
        l = np.array(res[0])
        u = np.array(res[1])
        a = u.dot(l)
        if iteration == 0:
            for i in range(length):
                previous[i] = a[i][i]
            iteration += 1
            continue
        for j in range(length):
            current[j] = abs(a[j][j] - previous[j])
            previous[j] = a[j][j]
        cur = max(current)
        iteration += 1
    return a, cur, iteration

# Записываем результат функции
# Верхнетреугольная матрица
# print('Получившаяся матрица: ')
# print(matrix2)


def find_eig_value(mat):
    eig_value = [0 for i in range(length)]
    for k in range(length):
        eig_value[k] = mat[k][k]
    eig_value = eig_value[::-1]
    return eig_value


# eig_value2 = find_eig_value()
# print('Получившиеся собственные числа: ')
# print(eig_value2)


# Создаем систему из (n-1)-ого уравнений, закладывая x(n) = 1


def character_polynomial(mat, eig_v):
    a = [[[0 for i in range(length - 1)] for j in range(length - 1)] for k in range(length)]
    b = [[0 for i in range(length - 1)] for k in range(length)]
    for k in range(length):
        for i in range(length - 1):
            b[k][i] = -mat[i][length - 1]
            for j in range(length - 1):
                a[k][i][j] = mat[i][j]
            a[k][i][i] = mat[i][i] - eig_v[k]
    return a, b


# Записываем результат характеристического многочлена
# Массив матриц для каждого собственного числа
# Массив свободных членов для каждого собственного числа


# Находим собственные вектора для каждого собственного числа


def find_eigenvectors(m_a, m_b):
    matrix_x = [[0 for i in range(length - 1)] for j in range(length)]  # Массив собственных векторов
    for k in range(length):
        matrix_x[k] = np.linalg.solve(m_a[k], m_b[k]).tolist()
        matrix_x[k].append(1)
        matrix_x[k] = matrix_x[k] / np.linalg.norm(matrix_x[k])
    return matrix_x


# # Собственные вектора для соответсвующих собственных чисел
# print('Получившиеся собственные вектора: ')
# print(eig_vector2)


# Находим значения синуса между заданным и получившимся собственнымыми векторами


def sin_phi(eig_v1, eig_v2):
    sin_phi = [0 for i in range(length)]
    for i in range(length):
        phi = np.arccos((eig_v1[i].dot(eig_v2[i])))
        sin_phi[i] = abs(np.sin(phi))
    sin_phi_ = max(sin_phi)
    return sin_phi_


def chart1_2_3():
    coefficients = [0 for i in range(length3)]
    sinus = [0 for i in range(length3)]
    rel = [0 for i in range(length3)]
    iteration_max = [0 for i in range(length3)]
    eig_vector1 = np.eye(length)
    for j in range(length3):
        coefficients[j] = 1 / (10 ** (j+1))
        mat_cf = do_matrix(coefficients[j])
        graph = lu_algorithm(mat_cf, 1e-14, 3000)
        eig_val2 = find_eig_value(graph[0])
        res1 = character_polynomial(mat_cf, eig_val2)
        eig_vector2 = find_eigenvectors(res1[0], res1[1])
        rel[j] = graph[1]
        iteration_max[j] = graph[2]
        sinus[j] = sin_phi(eig_vector1, eig_vector2)
    return coefficients, sinus, rel, iteration_max


result2 = chart1_2_3()
coefficient = result2[0]
# print(coefficient)
sinus_vectors = result2[1][::-1]
# print(sinus_vectors)
relative_error = result2[2]
# print(relative_error)
iteration_max_ = result2[3]
# print(iteration_max_)



plt.figure(1)
plt.grid()
plt.xscale('log')
plt.xlabel('Коэффициент отделимости собсвтенных чисел')
plt.ylabel('Синус угла между собственными векторами')
plt.title("График зависимости синуса угла между собственными векторами от отделимости собственных чисел")
plt.plot(coefficient, sinus_vectors)

plt.figure(2)
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Коэффициент отделимости собсвтенных чисел')
plt.ylabel('Относительная погрешность')
plt.title("График зависимости относительной погрешности от отделимости собственных чисел")
plt.plot(coefficient, relative_error)

plt.figure(3)
plt.grid()
plt.xscale('log')
plt.xlabel('Коэффициент отделимости собсвтенных чисел')
plt.ylabel('Максимальное количество итераций')
plt.title("График зависимости максимального количества итераций от отделимости собственных чисел")
plt.plot(coefficient, iteration_max_)


def chart4(p):
    rel = [0 for i in range(100)]
    iteration_max = [0 for i in range(100)]
    for j in range(100):
        iteration_max[j] = j + 101
        mat_cf = do_matrix(p)
        graph = lu_algorithm(mat_cf, 1e-16, iteration_max[j])
        rel[j] = graph[1]

    return rel, iteration_max


result3_1 = chart4(1)
relative1 = result3_1[0]
number_iter1 = result3_1[1]

result3_2 = chart4(1e-3)
relative2 = result3_2[0]
number_iter2 = result3_2[1]

plt.figure(4)
plt.grid()
plt.yscale('log')
plt.xlabel('Номер итерации')
plt.ylabel('Относительная погрешность')
plt.title("График зависимости относительной погрешности от номера итерации")
plt.plot(number_iter1, relative1, label='Отделимость собсвтенных чисел = 1')
plt.plot(number_iter2, relative2, label='Отделимость собсвтенных чисел = 1e-3')
plt.legend()


def chart5(p):
    rel = [0 for i in range(10)]
    iteration_max = [0 for i in range(10)]
    for j in range(10):
        rel[j] = 1 / (10 ** (j + 3))
        mat_cf = do_matrix(p)
        graph = lu_algorithm(mat_cf, rel[j], 300)
        iteration_max[j] = graph[2]

    return rel, iteration_max


result4_1 = chart5(1)
rel_err1 = result4_1[0]
iter_rel1 = result4_1[1]

result4_2 = chart5(1e-3)
rel_err2 = result4_2[0]
iter_rel2 = result4_2[1]

plt.figure(5)
plt.grid()
plt.xscale('log')
plt.xlabel('Задаваемая точность')
plt.ylabel('Максимальное количество итераций')
plt.title("График зависимости максимального количества итераций от задаваемой точности")
plt.plot(rel_err1, iter_rel1, label='Отделимость собсвтенных чисел = 1')
plt.plot(rel_err2, iter_rel2, label='Отделимость собсвтенных чисел = 1e-3')
plt.legend()


def chart6():
    rel = [0 for i in range(length4)]
    for j in range(length4):
        mat_cf = do_matrix(1)
        mat_dm = np.array(deforming_matrix(mat_cf))
        graph = lu_algorithm(mat_dm[j], 1e-16, 200)
        rel[j] = 1 - (lu_algorithm(mat_dm[0], 1e-16, 200)[1] / graph[1])
    return rel


error = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

plt.figure(6)
plt.grid()
plt.xlabel('Изменение исходной матрицы в процентах')
plt.ylabel('Относительная погрешность относительно исходной матрицы')
plt.title("График зависимости относительной погрешности от изменения исходной матрицы")
plt.plot(error, chart6(), label='Отделимость собсвтенных чисел = 1')
plt.legend()
plt.show()