""" Plotting routine for StokvecMatrix """

import pyarts3 as pyarts
import numpy as np
from . import Matrix

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.StokvecMatrix,
    *,
    fig=None,
    ax=None,
    xgrid: pyarts.arts.Vector | None = None,
    ygrid: pyarts.arts.Vector | None = None,
    component: pyarts.arts.Stokvec = pyarts.arts.Stokvec("I"),
    **kwargs
):
    """Plot a Stokes vector matrix as a 2D heatmap for a specific component.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a 2D radiance field
        nx, ny = 40, 30
        x = np.linspace(0, 2*np.pi, nx)
        y = np.linspace(0, 2*np.pi, ny)
        X, Y = np.meshgrid(x, y, indexing='ij')

        stokes_mat = pyarts.arts.StokvecMatrix(np.zeros((nx, ny, 4)))
        stokes_mat[:, :, 0] = np.sin(X) * np.cos(Y)  # I component

        pyarts.plots.StokvecMatrix.plot(stokes_mat)

    Parameters
    ----------
    data : ~pyarts3.arts.StokvecMatrix
        A matrix of Stokes vectors (4 components each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : ~pyarts3.arts.Vector | None = None,
        X-axis values. If None, uses column indices. Defaults to None.
    ygrid : ~pyarts3.arts.Vector | None = None,
        Y-axis values. If None, uses row indices. Defaults to None.
    component : ~pyarts3.arts.Stokvec, optional
        Which Stokes component to plot. Defaults to I.
    **kwargs
        Additional keyword arguments passed to imshow()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """

    return Matrix.plot(np.einsum("ijk,k->ij", data, component),
                       fig=fig,
                       ax=ax,
                       xgrid=xgrid,
                       ygrid=ygrid,
                       **kwargs)
