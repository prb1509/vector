# Vector      

Vector is a small Python package for basic vector operations. It works not only for vectors in a cartesian coordinate system, but also for polar, spherical, and cylindrical coordinates.

## Installation

    pip install git+https://github.com/prb1509/vector

## Some Notes on Implementation and Notation

- The curvilinear coordinates all follow ISO standard. 
For example, the angular coordinate is given by $\phi$ in polar and cylindrical coordinates, azimuthal angle is denoted $\phi$ and polar angle is denoted $\theta$ in spherical coordinates.

- In addition, the parameters for curvilinear vectors do not correspond to the basis vectors of the system. Consider the following snippet:

```python
    s = spherical(5,1,-1)
```

- This is not the same as $s = 5 \hat{r} + \hat{\phi} - \hat{\theta}$.
- Rather these are the coordinates corresponding to a radial vector of length 5, rotated by 1 radian in the xy plane, and rotated -1 radian along the z axis.

## License

[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)