import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d
np.set_printoptions(precision=3, linewidth=500, suppress=True)

def interpolate_frame(frame, resolution=(100,100), kind='linear'):
    dim_x, dim_y = frame.shape
    if np.iscomplexobj(frame):
        re_frame = interpolate_frame(np.real(frame),
                                     resolution=resolution,
                                     kind=kind)
        im_frame = interpolate_frame(np.imag(frame),
                                     resolution=resolution,
                                     kind=kind)
        return re_frame + 1j*im_frame
    else:
        f = interp2d(np.arange(dim_x), np.arange(dim_y), np.nan_to_num(frame), kind=kind)
        frame = f(np.linspace(0, dim_x-1, resolution[0]),
                  np.linspace(0, dim_y-1, resolution[1]))
    return frame

if __name__ == '__main__':
    A = np.array([[np.nan, np.nan, np.nan, 1., 3.],
                  [np.nan, np.nan, 1., 3., 3.],
                  [np.nan, 2., 1., 2., 4.],
                  [1., 3., 2., 0., 2.],
                  [1., 1., 3., 4., 5.]])

    B = copy.deepcopy(A)
    C = copy.deepcopy(B)
    D = copy.deepcopy(B)
    D[np.isnan(D)==False] = 1.
    D[np.isnan(D)] = 0.

    for i, row in enumerate(C):
        for j, item in enumerate(row):
            if np.isnan(item):
                a = B[i, j+1]  if j < len(row)-1 else np.nan
                b = B[i+1, j]  if i < len(C)-1 else np.nan
                c = B[i, j-1]  if j > 0 else np.nan
                d = B[i-1, j]  if i > 0 else np.nan
                arr = np.array([a, b, c, d])
                if len(arr) > 0:
                    C[i, j] = np.nanmean(arr)

    mask_cri = 0.5
    C_interpol = interpolate_frame(C)
    D_interpol = interpolate_frame(D)
    C_interpol[D_interpol<mask_cri] = np.nan
    D_interpol[D_interpol<mask_cri] = 0.

    print(A)
    print(B)
    print(C)
    print(D)
    print(C_interpol)
    print(D_interpol)

    fig, axs = plt.subplots(2,2)
    axs[0, 0].imshow(B, vmin=0.0, vmax=5.0)
    axs[0, 1].imshow(C_interpol, vmin=0.0, vmax=5.0)
    axs[1, 0].imshow(D, vmin=0.0, vmax=5.0)
    axs[1, 1].imshow(D_interpol, vmin=0.0, vmax=5.0)
    plt.show()
