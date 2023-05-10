# Author(s): Compas/Compas FEA Team, Marius  Weber (ETHZ, HSLU T&A)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from math import log



__all__ = [
    'Material',
    'Concrete',
    'ConcreteSmearedCrack',
    'ConcreteDamagedPlasticity',
    'ElasticIsotropic',
    'Stiff',
    'MPCStiff',    
    'ElasticOrthotropic',
    'ElasticPlastic',
    # 'ThermalMaterial',
    'Steel',
    'CMMUsermat'
]


class Material(object):
    """Initialises base Material object.

    Parameters
    ----------
    name : str
        Name of the Material object.

    Attributes
    ----------
    name : str
        Name of the Material object.

    """

    def __init__(self, name):

        self.__name__ = 'Material'
        self.name = name
        self.attr_list = ['name']

    def __str__(self):
        print('\n')
        print('compas_fea {0} object'.format(self.__name__))
        print('-' * (len(self.__name__) + 18))

        for attr in self.attr_list:
            print('{0:<11} : {1}'.format(attr, getattr(self, attr)))

        return ''

    def __repr__(self):
        return '{0}({1})'.format(self.__name__, self.name)


# ==============================================================================
# linear elastic
# ==============================================================================

class ElasticIsotropic(Material):
    """Elastic, isotropic and homogeneous material.

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].
    tension : bool
        Can take tension.
    compression : bool
        Can take compression.

    """

    def __init__(self, name, E, v, p, tension=True, compression=True):
        Material.__init__(self, name=name)

        self.__name__ = 'ElasticIsotropic'
        self.name = name
        self.E = {'E': E}
        self.v = {'v': v}
        self.G = {'G': 0.5 * E / (1 + v)}
        self.p = p
        self.tension = tension
        self.compression = compression
        self.attr_list.extend(['E', 'v', 'G', 'p', 'tension', 'compression'])


class Stiff(ElasticIsotropic):
    """Elastic, very stiff and massless material.

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].

    """

    def __init__(self, name, E=10**13):
        ElasticIsotropic.__init__(self, name=name, E=E, v=0.3, p=10**(-1))

        self.__name__ = 'Stiff'

class MPCStiff(ElasticIsotropic):
    """Elastic, very stiff and massless material for MPC

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].

    """

    def __init__(self, name, E=10**13):
        ElasticIsotropic.__init__(self, name=name, E=E, v=0.3, p=10**(-1))

        self.__name__ = 'MPCStiff'        


class ElasticOrthotropic(Material):
    """Elastic, orthotropic and homogeneous material.

    Parameters
    ----------
    name : str
        Material name.
    Ex : float
        Young's modulus Ex in x direction [Pa].
    Ey : float
        Young's modulus Ey in y direction [Pa].
    Ez : float
        Young's modulus Ez in z direction [Pa].
    vxy : float
        Poisson's ratio vxy in x-y directions [-].
    vyz : float
        Poisson's ratio vyz in y-z directions [-].
    vzx : float
        Poisson's ratio vzx in z-x directions [-].
    Gxy : float
        Shear modulus Gxy in x-y directions [Pa].
    Gyz : float
        Shear modulus Gyz in y-z directions [Pa].
    Gzx : float
        Shear modulus Gzx in z-x directions [Pa].
    p : float
        Density [kg/m3].
    tension : bool
        Can take tension.
    compression : bool
        Can take compression.

    Notes
    -----
    - Can be created but is currently not implemented.

    """

    def __init__(self, name, Ex, Ey, Ez, vxy, vyz, vzx, Gxy, Gyz, Gzx, p, tension=True, compression=True):
        Material.__init__(self, name=name)

        self.__name__ = 'ElasticOrthotropic'
        self.name = name
        self.E = {'Ex': Ex, 'Ey': Ey, 'Ez': Ez}
        self.v = {'vxy': vxy, 'vyz': vyz, 'vzx': vzx}
        self.G = {'Gxy': Gxy, 'Gyz': Gyz, 'Gzx': Gzx}
        self.p = p
        self.tension = tension
        self.compression = compression
        self.attr_list.extend(['E', 'v', 'G', 'p', 'tension', 'compression'])


# ==============================================================================
# non-linear general
# ==============================================================================

class ElasticPlastic(Material):
    """Elastic and plastic, isotropic and homogeneous material.

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].
    f : list
        Plastic stress data (positive tension values) [Pa].
    e : list
        Plastic strain data (positive tension values) [-].

    Notes
    -----
    - Plastic stress--strain pairs applies to both compression and tension.

    """

    def __init__(self, name, E, v, p, f, e):
        Material.__init__(self, name=name)

        fc = [-i for i in f]
        ec = [-i for i in e]

        self.__name__ = 'ElasticPlastic'
        self.name = name
        self.E = {'E': E}
        self.v = {'v': v}
        self.G = {'G': 0.5 * E / (1 + v)}
        self.p = p
        self.tension = {'f': f, 'e': e}
        self.compression = {'f': fc, 'e': ec}
        self.attr_list.extend(['E', 'v', 'G', 'p', 'tension', 'compression'])


# ==============================================================================
# non-linear metal
# ==============================================================================

class Steel(Material):
    """Bi-linear steel with given yield stress.

    Parameters
    ----------
    name : str
        Material name.
    fy : float
        Yield stress [MPa].
    fu : float
        Ultimate stress [MPa].
    eu : float
        Ultimate strain [%].
    E : float
        Young's modulus E [GPa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].

    """

    def __init__(self, name, fy=355, fu=None, eu=20, E=210, v=0.3, p=7850):
        Material.__init__(self, name=name)

        E *= 10.**9
        fy *= 10.**6
        eu *= 0.01

        if not fu:
            fu = fy
        else:
            fu *= 10.**6

        ep = eu - fy / E
        f = [fy, fu]
        e = [0, ep]
        fc = [-i for i in f]
        ec = [-i for i in e]

        self.__name__ = 'Steel'
        self.name = name
        self.fy = fy
        self.fu = fu
        self.eu = eu
        self.ep = ep
        self.E = {'E': E}
        self.v = {'v': v}
        self.G = {'G': 0.5 * E / (1 + v)}
        self.p = p
        self.tension = {'f': f, 'e': e}
        self.compression = {'f': fc, 'e': ec}
        self.attr_list.extend(['fy', 'fu', 'eu', 'ep', 'E', 'v', 'G', 'p', 'tension', 'compression'])


# ==============================================================================
# non-linear timber
# ==============================================================================


# ==============================================================================
# non-linear masonry
# ==============================================================================


# ==============================================================================
# non-linear concrete
# ==============================================================================

class Concrete(Material):
    """Elastic and plastic-cracking Eurocode based concrete material.

    Parameters
    ----------
    name : str
        Material name.
    fck : float
        Characteristic (5%) 28 day cylinder strength [MPa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].
    fr : list
        Failure ratios.

    Notes
    -----
    - The concrete model is based on Eurocode 2 up to fck=90 MPa.

    """

    def __init__(self, name, fck, v=0.2, p=2400, fr=None):
        Material.__init__(self, name=name)

        de = 0.0001
        fcm = fck + 8
        Ecm = 22 * 10**3 * (fcm / 10.)**0.3
        ec1 = min(0.7 * fcm**0.31, 2.8) * 0.001
        ecu1 = 0.0035 if fck < 50 else (2.8 + 27 * ((98 - fcm) / 100.)**4) * 0.001

        k = 1.05 * Ecm * ec1 / fcm
        e = [i * de for i in range(int(ecu1 / de) + 1)]
        ec = [ei - e[1] for ei in e[1:]]
        fctm = 0.3 * fck**(2. / 3.) if fck <= 50 else 2.12 * log(1 + fcm / 10.)
        f = [10**6 * fcm * (k * (ei / ec1) - (ei / ec1)**2) / (1. + (k - 2) * (ei / ec1)) for ei in e]

        E = f[1] / e[1]
        ft = [1., 0.]
        et = [0., 0.001]

        if not fr:
            fr = [1.16, fctm / fcm]

        self.__name__ = 'Concrete'
        self.name = name
        self.fck = fck * 10.**6
        self.E = {'E': E}
        self.v = {'v': v}
        self.G = {'G': 0.5 * E / (1 + v)}
        self.p = p
        self.tension = {'f': ft, 'e': et}
        self.compression = {'f': f[1:], 'e': ec}
        self.fratios = fr
        self.attr_list.extend(['fck', 'fratios', 'E', 'v', 'G', 'p', 'tension', 'compression'])


class ConcreteSmearedCrack(Material):
    """Elastic and plastic, cracking concrete material.

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].
    fc : list
        Plastic stress data in compression [Pa].
    ec : list
        Plastic strain data in compression [-].
    ft : list
        Plastic stress data in tension [-].
    et : list
        Plastic strain data in tension [-].
    fr : list
        Failure ratios.

    """

    def __init__(self, name, E, v, p, fc, ec, ft, et, fr=[1.16, 0.0836]):
        Material.__init__(self, name=name)

        self.__name__ = 'ConcreteSmearedCrack'
        self.name = name
        self.E = {'E': E}
        self.v = {'v': v}
        self.G = {'G': 0.5 * E / (1 + v)}
        self.p = p
        self.tension = {'f': ft, 'e': et}
        self.compression = {'f': fc, 'e': ec}
        self.fratios = fr
        self.attr_list.extend(['E', 'v', 'G', 'p', 'tension', 'compression', 'fratios'])


class ConcreteDamagedPlasticity(Material):
    """Damaged plasticity isotropic and homogeneous material.

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].
    damage : list
        Damage parameters.
    hardening : list
        Compression hardening parameters.
    stiffening : list
        Tension stiffening parameters.

    """

    def __init__(self, name, E, v, p, damage, hardening, stiffening):
        Material.__init__(self, name=name)

        self.__name__ = 'ConcreteDamagedPlasticity'
        self.name = name
        self.E = {'E': E}
        self.v = {'v': v}
        self.G = {'G': 0.5 * E / (1 + v)}
        self.p = p
        self.damage = damage
        self.hardening = hardening
        self.stiffening = stiffening
        self.attr_list.extend(['E', 'v', 'G', 'p', 'damage', 'hardening', 'stiffening'])


# ==============================================================================
# thermal
# ==============================================================================

class ThermalMaterial(Material):
    """Class for thermal material properties.

    Parameters
    ----------
    name : str
        Material name.
    conductivity : list
        Pairs of conductivity and temperature values.
    p : list
        Pairs of density and temperature values.
    sheat : list
        Pairs of specific heat and temperature values.

    """

    def __init__(self, name, conductivity, p, sheat):
        Material.__init__(self, name=name)

        self.__name__ = 'ThermalMaterial'
        self.name = name
        self.conductivity = conductivity
        self.p = p
        self.sheat = sheat
        self.attr_list.extend(['p', 'conductivity', 'sheat'])

# ==============================================================================
# convetional CMM Usermat
# ==============================================================================

class CMMUsermat(Material):
    """Elastic, isotropic and homogeneous material.

    Parameters
    ----------
    name : str
        Material name.
    E : float
        Young's modulus E [Pa].
    v : float
        Poisson's ratio v [-].
    p : float
        Density [kg/m3].
    tension : bool
        Can take tension.
    compression : bool
        Can take compression.

    """

    def __init__(self, name, geo, concrete, reinf_1L, reinf_2L, reinf_3L, reinf_4L):
        Material.__init__(self, name=name)

        self.__name__ = 'CMMUsermat'
        self.name = name

        # Delete nur als workflow test jetzt noch drin
        self.E = {'E': 200000}
        self.v = {'v': 0.1}
        self.G = {'G': 0.5 * 200000 / (1 + 0.1)}
        self.p = 0.0000025

        # Schlussendlich mit diesen werten
        self.R_Rohr = {'R_Rohr': geo['R_Rohr'] }
        self.rho = {'rho': geo['rho'] }    
        self.oo = {'oo': geo['oo'] }    
        self.uu = {'uu': geo['uu'] }    
        
        self.beton = {'beton': concrete['beton'] }    
        self.fcc = {'fcc': concrete['fcc'] }    
        self.vc = {'vc': concrete['vc'] }    
        self.ecu = {'ecu': concrete['ecu'] }    
        self.k_E = {'k_E': concrete['k_E'] }    
        self.theta_b0 = {'theta_b0': concrete['theta_b0'] }    
        self.theta_b1 = {'theta_b1': concrete['theta_b1'] }    
        self.k_riss = {'k_riss': concrete['k_riss'] }    
        self.Entfestigung = {'Entfestigung': concrete['Entfestigung'] }    
        self.lambdaTS = {'lambdaTS': concrete['lambdaTS'] }    
        self.srmx = {'srmx': concrete['srmx'] }    
        self.srmy = {'srmy': concrete['srmy'] }    
        self.Begrenzung = {'Begrenzung': concrete['Begrenzung'] }    
        self.KritQ = {'KritQ': concrete['KritQ'] }    
        self.winkelD = {'winkelD': concrete['winkelD'] }    
        self.k_vr = {'k_vr': concrete['k_vr'] }    
        self.fswy = {'fswy': concrete['fswy'] }   

        self.stahl1 = {'stahl1': reinf_1L['stahl'] }    
        self.zm1 = {'zm1': reinf_1L['zm'] }    
        self.fsy1 = {'fsy1': reinf_1L['fsy'] }    
        self.fsu1 = {'fsu1': reinf_1L['fsu'] }    
        self.esu1 = {'esu1': reinf_1L['esu'] }    
        self.esv1 = {'esv1': reinf_1L['esv'] }    
        self.Es1 = {'Es1': reinf_1L['Es'] }    
        self.ka1 = {'ka1': reinf_1L['ka'] }    
        self.kb1 = {'kb1': reinf_1L['kb'] }    
        self.kc1 = {'kc1': reinf_1L['kc'] }    
        self.as1 = {'as1': reinf_1L['as'] }    
        self.dm1 = {'dm1': reinf_1L['dm'] }    
        self.psi1 = {'psi1': reinf_1L['psi'] }    

        self.stahl2 = {'stahl2': reinf_2L['stahl'] }    
        self.zm2 = {'zm2': reinf_2L['zm'] }    
        self.fsy2 = {'fsy2': reinf_2L['fsy'] }    
        self.fsu2 = {'fsu2': reinf_2L['fsu'] }    
        self.esu2 = {'esu2': reinf_2L['esu'] }    
        self.esv2 = {'esv2': reinf_2L['esv'] }    
        self.Es2 = {'Es2': reinf_2L['Es'] }    
        self.ka2 = {'ka2': reinf_2L['ka'] }    
        self.kb2 = {'kb2': reinf_2L['kb'] }    
        self.kc2 = {'kc2': reinf_2L['kc'] }    
        self.as2 = {'as2': reinf_2L['as'] }    
        self.dm2 = {'dm2': reinf_2L['dm'] }    
        self.psi2 = {'psi2': reinf_2L['psi'] }    

        self.stahl3 = {'stahl3': reinf_3L['stahl'] }    
        self.zm3 = {'zm3': reinf_3L['zm'] }    
        self.fsy3 = {'fsy3': reinf_3L['fsy'] }    
        self.fsu3 = {'fsu3': reinf_3L['fsu'] }    
        self.esu3 = {'esu3': reinf_3L['esu'] }    
        self.esv3 = {'esv3': reinf_3L['esv'] }    
        self.Es3 = {'Es3': reinf_3L['Es'] }    
        self.ka3 = {'ka3': reinf_3L['ka'] }    
        self.kb3 = {'kb3': reinf_3L['kb'] }    
        self.kc3 = {'kc3': reinf_3L['kc'] }    
        self.as3 = {'as3': reinf_3L['as'] }    
        self.dm3 = {'dm3': reinf_3L['dm'] }    
        self.psi3 = {'psi3': reinf_3L['psi'] }   

        self.stahl4 = {'stahl4': reinf_4L['stahl'] }    
        self.zm4 = {'zm4': reinf_4L['zm'] }    
        self.fsy4 = {'fsy4': reinf_4L['fsy'] }    
        self.fsu4 = {'fsu4': reinf_4L['fsu'] }    
        self.esu4 = {'esu4': reinf_4L['esu'] }    
        self.esv4 = {'esv4': reinf_4L['esv'] }    
        self.Es4 = {'Es4': reinf_4L['Es'] }    
        self.ka4 = {'ka4': reinf_4L['ka'] }    
        self.kb4 = {'kb4': reinf_4L['kb'] }    
        self.kc4 = {'kc4': reinf_4L['kc'] }    
        self.as4 = {'as4': reinf_4L['as'] }    
        self.dm4 = {'dm4': reinf_4L['dm'] }    
        self.psi4 = {'psi4': reinf_4L['psi'] }                      
        
        self.attr_list.extend(['E', 'v', 'G', 'p', 'R_Rohr', 'rho', 'oo', 'uu', 'beton', 'fcc', 'vc', 'ecu', 'k_E', 'theta_b0', 'theta_b1', 'k_riss', 'Entfestigung', 'lambdaTS', 'srmx', 'srmy', 'Begrenzung', 'KritQ', 'winkelD', 'k_vr', 'fswy', 'stahl1', 'zm1', 'fsy1', 'fsu1', 'esu1', 'esv1', 'Es1', 'ka1', 'kb1', 'kc1', 'as1', 'dm1', 'psi1', 'stahl2', 'zm2', 'fsy2', 'fsu2', 'esu2', 'esv2', 'Es2', 'ka2', 'kb2', 'kc2', 'as2', 'dm2', 'psi2', 'stahl3', 'zm3', 'fsy3', 'fsu3', 'esu3', 'esv3', 'Es3', 'ka3', 'kb3', 'kc3', 'as3', 'dm3', 'psi3', 'stahl4', 'zm4', 'fsy4', 'fsu4', 'esu4', 'esv4', 'Es4', 'ka4', 'kb4', 'kc4', 'as4', 'dm4', 'psi4'])        