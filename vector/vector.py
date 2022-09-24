import cmath
import math

class VectorError(Exception):
    pass


class vector_base:
    """A basic vector class. Parent class for specific vector 
    classes cartesian, polar etc to inherit from.
    """
    #define magic methods
    def __init__(self,*components):
        self.vector = components
        
        if type(self) in [cylindrical, spherical ] and len(components) !=3:
            raise VectorError("Cylindrical/Spherical vectors must have exactly 3 components ")
        elif type(self) == polar and len(components) != 2:
            raise VectorError("Polar vectors must have exactly 2 components ")

        if len(components) == 0:
            self.vector = ()

    def __len__(self):
        """Return len(self.vector)"""
        return len(self.vector)

    def __repr__(self):
       """ Return repr(self)"""
       return str(self.vector)
    
    def __getitem__(self,key):
        """ Return self.vector[key]"""
        return self.vector[key]

    def __ne__(self,other):
        """Return self != other"""
        if len(self) != len(other):
            return True
        else:
            for i,j in zip(self,other):
                if i!= j:
                    return True
            else:
                return False

    #define binary operations
    def __add__(self,other):
        """Return self + other"""
        #we want to only add vectors of the same class or to the zero vector class
        #and the vectors must have the same length
        if (len(self) == len(other)):    
            if type(self) == type(other):                
                if type(self) == cartesian:
                        return cartesian(*[self.vector[i] + other.vector[i]
                                     for i in range(len(self))])
           
                elif type(self) == polar:
                        return cartesian.to_polar(
                            polar.to_rect(self) + polar.to_rect(other))
                    
                elif type(self) == cylindrical:
                    return cartesian.to_cylindrical(
                        cylindrical.to_rect(self) + cylindrical.to_rect(other))

                elif type(self) == spherical:
                    return cartesian.to_spherical(
                        spherical.to_rect(self) + spherical.to_rect(other))

                elif type(self) == zero:
                    return self

            #cases when one vector is zero and other is not, but is still a vector
            elif type(self) != type(other):
                if type(self) == zero:
                    return other

                elif type(other) == zero:
                    return self

                else:
                    raise VectorError("Incompatible vector types. Both vectors must be of the same class.")
        else:
            raise VectorError("Vectors must have the same length.")
        

    def __sub__(self,other):
        """Return self - other"""
        return self + -other
    

    def __mul__(self,other):
        """Return self * other"""
        #we want to multiply only vectors of the same class, unless again one of them is the zero vector, and also 
        #of the same length
        #in addition, we also need to define scalar multiplication 
        if type(self) in (zero,cartesian,polar,cylindrical,spherical) and type(other) in (zero,cartesian,polar,cylindrical,spherical) and type(self) == type(other):
            if len(self) == len(other):
                if type(self) == cartesian:
                    return sum(self[i].conjugate() * other[i]
                       for i in range(len(self)))

                elif type(self) == polar:
                    return polar.to_rect(self) * polar.to_rect(other)

                elif type(self) == cylindrical:
                    return cylindrical.to_rect(self) * cylindrical.to_rect(other)

                elif type(self) == spherical:
                    return spherical.to_rect(self) * spherical.to_rect(other)

                elif type(self) == zero:
                    return 0
            
            else:
                raise VectorError("Vectors must have the same length.")

        elif type(self) in (zero,cartesian,polar,cylindrical,spherical) and type(other) in (zero,cartesian,polar,cylindrical,spherical) and type(self) != type(other):
            if len(self) == len(other):
                if type(self) == zero or type(other) == zero:
                    return 0

                else:
                    raise VectorError("Incompatible vector types. Both vectors must be of the same class.")

            else:
                raise VectorError("Vectors must have the same length.")

        else:
            if type(self) == cartesian:
                return cartesian(*[self[i] * other for i in range(len(self))])

            elif type(self) == polar:
                 return cartesian.to_polar(polar.to_rect(self) * other)

            elif type(self) == cylindrical:
                return cartesian.to_cylindrical(
                    cylindrical.to_rect(self) * other)

            elif type(self) == spherical:
                return cartesian.to_spherical(spherical.to_rect(self) * other)

            elif type(self) == zero:
                return 0
            
    def __rmul__(self,other):
        """return self * other""" 
        return self * other
            
    def __truediv__(self,other):
        """Return self / other"""
        if type(other) in (int,float,complex):
            if type(self) == cartesian:
                return cartesian(*[self.vector[i] / other
                              for i in range(len(self))])

            elif type(self) == polar:
                return cartesian.to_polar(polar.to_rect(self) * 1 / other)

            elif type(self) == cylindrical:
                return cartesian.to_cylindrical(
                    cylindrical.to_rect(self) * 1 / other)

            elif type(self) == spherical:
                return cartesian.to_spherical(
                    spherical.to_rect(self) * 1 / other)

            elif type(self) == zero:
                return self
            
        else:
            raise VectorError("Vector can only be divided by int,float,or complex")

    def __floordiv__(self,other):
        """Return self // other"""
        if self.notreal():
             raise TypeError("Can't take floor of non real vector")

        if type(other) == int or type(other) == float:
            if type(self) == cartesian and cartesian.notreal(self):
                return cartesian(*[self.vector[i] // other
                              for i in range(len(self))])

            elif type(self) == polar:
                return cartesian.to_polar(polar.to_rect(self) * 1 // other)

            elif type(self) == cylindrical:
                return cartesian.to_cylindrical(
                    cylindrical.to_rect(self) * 1 // other)

            elif type(self) == spherical:
                return cartesian.to_spherical(
                    spherical.to_rect(self) * 1 // other)

            elif type(self) == zero:
                return self

        else:
            raise TypeError("Floor division defined for int,float denominator.")

    #define unary methods
    def __neg__(self):
        """Return -self"""
        return -1*self

    def __pos__(self):
        """Return +self"""
        return self

    #define the methods which are not unique to any vector class
    def inner(v,u):
        """Return inner product for vectors in given coordinates.
        For any two vectors v,u in C^n, the inner product <v,u> is defined as
        the product of v* and u, where v* is the conjugate transpose of v.
        
        Parameters:
        v: cartesian, zero, polar, spherical, cylindrical
        The first vector.

        u: cartesian, zero, polar, spherical, cylindrical
        The second vector.
        """
        return v * u

    def isreal(v):
        """Check if vector has only real components

        Parameters:
        v: cartesian, zero, polar, spherical, cylindrical
        The vector which is being checked."""
        for component in v:
            if type(component) not in [int,float]:
                return False
        else:
            return True

    def notreal(v):
        """Boolean inverse of isreal

        Parameters:
        v: cartesian, zero, polar, spherical, cylindrical
        The vector which is being checked."""
        return not v.isreal()


class zero(vector_base):
    """ Vector with 0 for all components; ie, zero vector."""
    def __init__(self,dimension):
        """Return zero vector of given dimension.
        Parameters:
        dimension: int
        The dimension of desired zero vector
        """
        self.vector = vector_base(*[0 for i in range(dimension)])

    
class cartesian(vector_base):
    """Cartesian vector of any dimension."""
    def norm(v):
        """Return norm of vector.

        Parameters:
        v: cartesian
        The vector whose norm is to be found.
        """
        return (v * v).real**0.5   

    def normalize(v):
        """Return normalized unit vector.

        Parameters:
        v: cartesian
        The vector to be normalized.
        """
        v.vector = v / v.norm()
        return v

    def angle(v,u,radians = True):
        """Return angle between vectors.

        Parameters:
        v: cartesian
        First vector.

        u: cartesian
        Second vector.
        """
        if radians:
            return math.acos((v * u).real / (v.norm() * u.norm()))

        else:
            return math.degrees(math.acos((v * u).real / (v.norm() * u.norm())))

    def cross(v,u):
        """Return vector cross prodcut for 3d or 7d vectors.

        Parameters:
        v: cartesian
        First vector of cross product

        u: cartesian
        Second vector of cross product
        """
        if len(v) == len(u) == 3:
            return cartesian(*(
                [v[1] * u[2] - v[2] * u[1],
                v[2] * u[0] - v[0] * u[2],
                 v[0] * u[1] - v[1] * u[0]]))

        elif len(v) == len(u) == 7:
            return cartesian(*(
                [v[1] * u[3] - v[3] * u[1] + v[2] * u[6] - v[6] * u[2] + v[4] * u[5] - v[5] * u[4],
                v[2] * u[4] - v[4] * u[2] + v[3] * u[0] - v[0] * u[3] + v[5] * u[6] - v[6] * u[5],
                v[3] * u[5] - v[5] * u[3] + v[4] * u[1] - v[1] * u[4] + v[6] * u[0] - v[0] * u[6],
                v[4] * u[6] - v[6] * u[4] + v[5] * u[2] - v[2] * u[5] + v[0] * u[1] - v[1] * u[0],
                v[5] * u[0] - v[0] * u[5] + v[6] * u[3] - v[3] * u[6] + v[1] * u[2] - v[2] * u[1],
                v[6] * u[1] - v[1] * u[6] + v[0] * u[4] - v[4] * u[0] + v[2] * u[3] - v[3] * u[2],
                v[0] * u[2] - v[2] * u[0] + v[1] * u[5] - v[5] * u[1] + v[3] * u[4] - v[4] * u[3],]))

        else:
            raise VectorError("Cross product defined for 3 and 7 dimensional vectors only.") 

    def rotate_2d(v,angle, radians = False):
        """Return 2D rotation of vector.

        Parameters:
        v: cartesian
        The vector to be rotated.

        angle: int,float
        Angle of rotation

        radians: bool,optional
        If True, input angle unit will be considered as radian
        """
        if len(v) != 2:
            raise ValueError("2D rotaion is defined for 2D vectors only")

        else:
            #multiply vector by the rotation matrix
            #[[cos(x), -sin(x)],[sin(x),cos(x)]] for angle x
            if radians:
                return cartesian(*(v[0] * math.cos(angle)- v[1] * math.sin(angle),
                            v[0] * math.sin(angle) + v[1] * math.cos(angle)))
                            
            else:
                angle = math.radians(angle)
                return cartesian(*(v[0] * math.cos(angle)- v[1] * math.sin(angle),
                            v[0] * math.sin(angle) + v[1] * math.cos(angle)))

    def proj(v,u):
        """Return projection of v on u.

        Parameters:
        v: cartesian
        The vector projected.

        u: cartesian
        The vector to be projected on.
        """
        #projection defined as (v.u)u/|u|**2
        if u.norm() == 0:
            raise VectorError("Cannot project onto a zero vector.")

        else:
            return (v * u) * u / (u.norm()**2)

    def conjugate(v):
        """Return conjugate of each component of the vector.

        Parameters:
        v: cartesian
        The vector whose components' conjugate is to be returned
        """
        return cartesian(*[v[i].conjugate() for i in range(len(v))])

    def imag(v):
        """Return imaginary part of each component of the vector.

        Parameters:
        v: cartesian
        The vector whose components' imaginary part is to be returned
        """
        return cartesian(*[v[i].imag for i in range(len(v))])

    def real(v):
        """Return real part of each component of the vector.

        Parameters:
        v: cartesian
        The vector whose components' real part is to be returned
        """
        return cartesian(*[v[i].real for i in range(len(v))])

    def dist(v,u):
        """Return distance between cartesian vectors.

        Parameters:
        v: cartesian
        The first vector.

        u: cartesian
        The second vector.
        """
        return math.sqrt(sum(cartesian(*[(v[i]-u[i])**2 for i in range(len(v))])))

    
    def to_polar(v):
        """Convert cartesian vector to polar vector.

        Parameters:
        v: cartesian 
        The vector to be converted to polar."""
        if len(v) == 2:
            r = math.sqrt(v[0]**2 + v[1]**2)
            phi = math.atan2(v[1],v[0])
            return polar(r,phi)

        else:
            raise VectorError("Vector must be 2 dimensional.")


    def to_cylindrical(v):
        """Convert cartesian vector to cylindrical vector.

        Parameters:
        v: cartesian
        The vector to be converted to cylindrical."""
        if len(v) == 3:
            rho = math.sqrt(v[0]**2 + v[1]**2)
            phi = math.atan2(v[1],v[0])
            if phi < 0:
                phi += math.pi
            z = v[2]
            return cylindrical(rho,phi,z)

        else:
            raise VectorError("Vector must be 3 dimensional.")


    def to_spherical(v):
        """Convert cartesian vector to spherical vector.

        Parameters:
        v: cartesian
        The vector to be converted to spherical."""
        if len(v) == 3:
            r = math.sqrt(v[0]**2+v[1]**2 + v[2]**2)
            phi = math.atan2(v[1],v[0])
            if phi < 0:
                phi += math.pi
            theta = math.acos(v[2]/v.norm())
            return spherical(r,theta,phi)

        else:
            raise VectorError("Vector must be 3 dimensional.")
        

class polar(vector_base):
    """2d polar vector."""
    @property
    def r(self):
        return self[0]

    @property
    def phi(self):
        return self[1]
    
    def dist(v,u):
        """Return distance between polar vectors.

        Parameters:
        v: polar
        The first vector.

        u: polar
        The second vector.
        """
        return math.sqrt(v.r**2 + u.r**2 -
                         2 * v.r * u.r * math.cos(v.phi - u.phi))

    def to_rect(v):
        """Convert polar coordinate vector to cartesian vector.

        Parameters:
        v: polar
        The vector to be converted to rectangular coordinates."""
        return v.r * cartesian(math.cos(v.phi),math.sin(v.phi))


class spherical(vector_base):
    """3d spherical vector."""        
    @property
    def r(self):
        return self[0]

    @property
    #polar angle
    def theta(self):
        return self[1]

    @property
    #azimuthal angle
    def phi(self):
        return self[2]

    def dist(v,u):
        """Return distance between spherical vectors.

        Parameters:
        v: spherical
        The first vector.

        u: spherical
        The second vector.
        """
        return math.sqrt(v.r**2 + u.r**2 -
                         2 * v.r * u.r *
                         (math.sin(v.phi) * math.sin(u.phi) *
                          math.cos(v.theta - u.theta) + math.cos(v.phi) *
                         math.cos(u.phi)))
    
    def to_rect(v):
        """Convert spherical coordinate vector to cartesian vector.

        Parameters:
        v: spherical
        The vector to be converted to rectangular coordinates."""
        x = v.r * math.cos(v.theta) * math.sin(v.phi)
        y = v.r * math.sin(v.theta) * math.sin(v.phi)
        z = v.r * math.cos(v.phi)
        return cartesian(x,y,z)

    def to_cylindrical(v):
        """Convert spherical coordinate vector to cartesian vector.

        Parameters:
        v: spherical
        The vector to be converted to rectangular coordinates."""
        rho = v.r * math.sin(v.phi)
        z = v.r * math.cos(v.phi)
        return cylindrical(rho,v.theta,z)

    
class cylindrical(vector_base):
    """3d cylindrical vector."""                                
    @property
    def rho(self):
        return self[0]

    @property
    def phi(self):
        return self[1]

    @property
    def z(self):
        return self[2]

    def norm(v):
        """Return norm of cylindrical vector.

        Parameters:
        v: cylindrical
        The vector whose norm is to be found.
        """
        return math.sqrt(v.rho**2 + v.z**2)

    def dist(v,u):
        """Return distance between cylindrical vectors.

        Parameters:
        v: cylindrical
        The first vector.

        u: cylindrical
        The second vector.
        """
        return math.sqrt(v.rho**2 + u.rho**2 -
                         2 * v.rho * u.rho * math.cos(v.phi - u.phi) +
                         (v.z - u.z)**2)
    
    def to_rect(v):
        """Convert cylindrical coordinate vector to cartesian vector.

        Parameters:
        v: cylindrical
        The vector to be converted to rectangular coordinates."""
        x = v.rho * math.cos(v.phi)
        y = v.rho * math.sin(v.phi)
        return cartesian(x,y,v.z)

    def to_spherical(v):
        """Convert cylindrical coordinate vector to cartesian vector.

        Parameters:
        v: cylindrical
        The vector to be converted to rectangular coordinates."""
        r = math.sqrt(v.rho**2 + v.z**2)
        theta = math.acos(v.z/r)
        return spherical(r,theta,v.phi)