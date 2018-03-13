from scipy.optimize import curve_fit
import numpy as np

def fit_power_law(x,y, T_c):
    def fit_func(T, a, b):
        return a * (T-T_c)**b
    params = curve_fit(fit_func, x, y)
    return params[0];

#essential test
# import matplotlib.pyplot as plt
# x = np.linspace(0.201, 1, 1000);
# y = np.abs(x-0.2)**0.1;
# plt.plot(x,y);
# T_c = 0.2
# params = fit_power_law(x,y,T_c);
# print(params);
# reconstruct = params[0]*(x-T_c)**params[1]
# plt.plot(x, reconstruct)
# plt.show()


