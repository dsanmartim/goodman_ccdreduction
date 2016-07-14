import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import gridspec
from numpy.polynomial.chebyshev import chebfit, chebval
from numpy.polynomial.legendre import legfit, legval
from scipy.interpolate import interp1d


def fitting(x, y, function='Legendre', order=5):
    """
    :rtype : array

    data: array

    niter: int or 'None'
    low_rej = int
        A factor of the data.std()
    high_rej = int
        A factor of the data.std()
    """

    # prepare for masking arrays - 'conventional' arrays won't do it
    # data = np.ma.array(data)

    # Assuring that the output yfit has always same number of points of the input
    global yfit
    x = np.arange(0, x.size, 1)

    if function.lower() == 'legendre':
        yfit = chebval(x, (chebfit(x, y, deg=order)))
    elif function.lower() == 'chebychev':
        yfit = legval(x, (legfit(x, y, deg=order)))
    elif function.lower() == 'spline3':
        nsum = 10

        y_resampled = np.asarray([np.median(y[i:i + nsum]) for i in range(0, len(y) - len(y) % nsum, nsum)])
        x_resampled = np.linspace(0, y.size, y_resampled.size)

        # Masking invalid elements (e.g. with nan)
        y_resampled = ma.masked_invalid(y_resampled)

        # Fitting
        f = interp1d(x_resampled, y_resampled, kind=order, bounds_error=False, fill_value=0.0)

        # Function with the original x array
        yfit = np.asarray(f(x))

        # TODO eliminar esse plots quando a rotina estiver funcionando
        # plotting Spline to see results
        plt.plot(x, y, 'k-', label='Original')
        plt.plot(x_resampled, y_resampled, 'bo', label='Resampled')
        plt.plot(x, yfit, label='Fit')
        plt.legend(loc='best')
        plt.show()

    # Calculating chi2
    npar = order + 1
    sigma2 = 1. / (y.size - 1.) * ((y - y.mean()) ** 2).sum()
    chi2_dof = ((y - yfit) ** 2 / sigma2).sum() / (y.size - (npar))

    return yfit, chi2_dof


def sigma_clipping(data, low_rej=3.0, high_rej=3.0, niter=5):
    # lower = - low_rej * data.std()
    # upper = + high_rej * data.std()
    # Just the good data (it masks vector outside the interval -lower and upper)
    # good_data = np.ma.masked_outside(data, lower, upper)
    # Just the clipped data
    # bad_data = np.ma.masked_inside(data, lower, upper)

    good_data = sigma_clip(data, sigma_lower=low_rej, sigma_upper=high_rej, iters=1, cenfunc=ma.mean)
    bad_data = good_data.mask * data
    bad_data = ma.masked_where(bad_data == 0.0, bad_data)

    # Creating masks
    good_mask = good_data / good_data
    bad_mask = bad_data / bad_data

    return good_mask, bad_mask


# Fitting Sline3
def fit_spline3(y, x, order=3, nsum=3):
    y_resampled = [np.median(y[i:i + nsum]) for i in range(0, len(y) - len(y) % nsum, nsum)]
    x_resampled = np.linspace(0, len(y), len(y_resampled))

    # Fitting
    f = interp1d(x_resampled, y_resampled, kind=order, bounds_error=True)

    # Return function to be constructed with any other x array
    return f


# Local Minima and Maxima
def local_minmax(data, nmin=2, nmax=2):
    # Identifying indices of local minima-maxima points
    id_min = (np.gradient(np.sign(np.gradient(data))) > 0).nonzero()[0]  # index of local min
    id_max = (np.gradient(np.sign(np.gradient(data))) < 0).nonzero()[0]  # index of local max

    # Taking values at min/max points
    list_min, list_max = data[id_min], data[id_max]

    # Sorting minima-maxima values (bigger --> lower)
    list_min, id_min = (list(p) for p in zip(*sorted(zip(list_min, id_min), reverse=False)))
    list_max, id_max = (list(p) for p in zip(*sorted(zip(list_max, id_max), reverse=True)))

    # Taking the desired number of local minima-maxima points
    list_min, list_max, id_min, id_max = list_min[0:nmin], list_max[0:nmax], id_min[0:nmin], id_max[0:nmax]

    return list_min, list_max, id_min, id_max


def trim_slitedge(flat):
    # Getting input data
    ccddata = fits.getdata(flat, ignore_missing_end=True)

    # Collapse flat in the dispersion direction
    flat_collapsed = fits.getdata(flat, ignore_missing_end=True).sum(axis=1) / ccddata.shape[1]
    lines = np.arange(0, flat_collapsed.size, 1)

    # Excluding first pixels in the spatial direction
    cut = 3
    c_flat = flat_collapsed[cut:-cut]
    c_lines = np.arange(0, c_flat.size, 1)

    # Fittin cubic spline. It's working very well with order=5, nsum=2
    func_splin3 = fit_spline3(c_flat, c_lines, order=5, nsum=2)
    smooth_flat = func_splin3(c_lines)

    # Compute 1st and flat smoothed
    dy = np.gradient(smooth_flat)
    dy2 = np.gradient(dy)

    # Regions to compute local minina-maxima
    # Region one: it represent first 40 percent of all data
    # Region two: ... last 40%
    pixa, pixb = int(len(c_flat) * 0.4), int(len(c_flat) * 0.6)
    dy2_one, dy2_two = dy2[0:pixa], dy2[pixb:]

    # Reg. 1: Compute local min/max of the 2nd derivative
    list_min_1, list_max_1, id_min_1, id_max_1 = local_minmax(dy2_one, nmin=1, nmax=1)
    list_min_2, list_max_2, id_min_2, id_max_2 = local_minmax(dy2_two, nmin=1, nmax=1)

    # Indice have to be reshifted to the original indices of the function dy2
    id_min_2 = np.array(id_min_2) + pixb

    # Slit edges are the local maxima/minima 1/2 [accounting the cutted pixels]
    slit_1, slit_2 = int(np.array(id_min_1) + cut), int(np.array(id_min_2) + cut)

    print slit_1, slit_2

    import matplotlib.pyplot as plt
    c_lines += cut
    plt.plot(lines, flat_collapsed, 'k-', label='Flat Collapsed')
    plt.plot(lines[slit_1:slit_2], flat_collapsed[slit_1:slit_2], 'r-', label = 'Cutted Flat')
    plt.plot(c_lines, dy, 'g-', label="Dy/dx")
    plt.plot(c_lines, dy2, 'y-', label="Dy2/dx")
    plt.plot(slit_1, list_min_1, 'bo', label='Slit Edge 1 ')
    plt.plot(slit_2, list_min_2, 'ro', label='Slit Edge 2')
    plt.xlim(lines.min() - 50, lines.max() + 50)
    plt.legend(loc='best')
    plt.show()

    return slit_1, slit_2


def cor_flat(flat, function='Legendre', order=5, clipping='sigma_clipping', niter=1, low_rej=3.0, high_rej=3.0):
    """
    Remove arc features of flat images from Goodman internal illumination. The ouptut file will be the
    input flat corrected.

    Parameters
    ----------
    flat : str
           Flat file
    function : str
           'Legendre' or 'Chebyshev'.
    order : int
           Maximum. order is defined by the respective built-in function
    clipping: str
           sigma_clipping, None
    """

    # Input and trimming data
    ccddata = fits.getdata(flat, ignore_missing_end=True)
    # ccddata = stats.trimboth(ccddata, 0.2, axis=0)

    # Collapsing 2D image into a single vector (a big aperture of all lines)
    spec = ccddata.sum(axis=0) / ccddata.shape[0]
    pix = np.arange(0, spec.size, 1)

    # Fitting function to the spectrum of the flat before any clipping
    firstfit, chi2_dof = fitting(pix, spec, function=function, order=order)

    # Buffering quantities before clipping
    firstspec = spec
    firstresidual = spec - firstfit

    # TODO Apagar ao final
    plt.plot(pix, firstspec, 'k-', label='1st Spec - Before Clip')
    plt.plot(pix, firstfit, 'g-', label='1st Fit - Before Clip')

    # Recovering the original spectrum ans residual (data - fitting)
    fit = firstfit
    residual = firstresidual

    if clipping.lower() == 'sigma_clipping' and niter > 0:

        niter += 1
        for iters in np.arange(1, niter, 1):
            good_mask, bad_mask = sigma_clipping(residual, low_rej=low_rej, high_rej=high_rej, niter=0)

            # Filtering spec for good data
            spec = spec * good_mask

            # Fitting after  clipping iteration
            fit, chi2_dof = fitting(pix, spec, function=function, order=order)

            print fit - firstfit
            # TODO Continuar checagem do codigo a partir daqui...

            # print fit - firstfit
            residual = spec - fit

            # plt.plot(pix, spec, 'ro', label='Spec After Clip')
            plt.plot(pix, firstspec * bad_mask, 'gx', label='Clipped')
            plt.plot(pix, fit + fit.mean() * 0.01, 'b-', label='Fit After  Clip')
            plt.legend(loc='best')
            plt.show()


    # Plotting things
    fig = plt.figure()
    plt.clf()

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex=ax0)

    firstresidual = 100 * (firstresidual / firstspec)
    residual = 100 * (residual / spec)

    if clipping.lower() == 'sigma_clipping' and niter > 0:
        ax0.plot(pix, firstspec * bad_mask, 'rx', label='Clipped', lw=1)
        ax1.plot(pix, firstresidual * bad_mask, 'rx', alpha=0.5, lw=2)

        # ax1.plot(pix, bad_res, 'x', color='r', alpha=0.5, lw=2)

    ax0.plot(pix, firstspec, 'k-', label='Flat Spectrum', lw=2)
    ax0.plot(pix, fit, 'g-', label=str(function.title()) + ' - Order ' + str(order))
    ax1.plot(pix, residual, '-', color='grey', alpha=0.5, lw=2)

    ax0.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=12)

    ax0.set_title('Chi2 = ' + str(chi2_dof))

    ax0.set_ylabel("ADU's", fontsize=12)
    ax1.set_xlabel('Column [pixel]', fontsize=12)
    ax1.set_ylabel('residual [%]', fontsize=11)

    ax0.set_xlim(pix[0] - 150, pix[-1] + 150)
    ax1.set_ylim(residual.min() - 4 * residual.std(), residual.max() + 4 * residual.std())
    plt.show()


flat_name = '/home/davidsanm/PyCharmProjects/GoodmanDataReduction/2016-03-20/RED/master_flat_600.fits'
#cor_flat(flat_name, function='Spline3', order=9, clipping='sigma_clipping', niter=0, low_rej=1.5, high_rej=1.5)

trim_slitedge(flat_name)
