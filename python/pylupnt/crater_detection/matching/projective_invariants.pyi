from __future__ import annotations
import numpy as np
from numpy import linalg as LA
from pylupnt.crater_detection.common.conics import conic_center
from pylupnt.crater_detection.common.conics import conic_matrix
from pylupnt.crater_detection.common.conics import matrix_adjugate
from pylupnt.crater_detection.common.conics import scale_det
from pylupnt.crater_detection.matching.utils import enhanced_pattern_shifting
from pylupnt.crater_detection.matching.utils import is_clockwise
from pylupnt.crater_detection.matching.utils import is_colinear
from pylupnt.crater_detection.matching.utils import np_swap_columns
from pylupnt.crater_detection.matching.utils import shift_nd
__all__ = ['CoplanarInvariants', 'LA', 'PermutationInvariant', 'conic_center', 'conic_matrix', 'enhanced_pattern_shifting', 'is_clockwise', 'is_colinear', 'matrix_adjugate', 'np', 'np_swap_columns', 'scale_det', 'shift_nd']
class CoplanarInvariants:
    @classmethod
    def from_detection_conics(cls, A_craters, crater_triads = None):
        ...
    @classmethod
    def match_generator(cls, A_craters = None, x_pix = None, y_pix = None, a_pix = None, b_pix = None, psi_pix = None, batch_size = 1, convert_to_radians = True, max_iter = 10000, sort_ij = True):
        """
        Generator function that yields crater triad and its associated projective invariants [1]. Triads are formed
                using Enhanced Pattern Shifting method [2]. Input craters can either be parsed as  parameterized ellipses
                (x_pix, y_pix, a_pix, b_pix, psi_pix) or as matrix representation of conic.
        
                Parameters
                ----------
                A_craters : np.ndarray
                    Crater detections in conic representation.
                x_pix, y_pix : np.ndarray
                    Crater center positions in image plane
                a_pix, b_pix : np.ndarray
                    Crater ellipse axis parameters in image plane
                psi_pix : np.ndarray
                    Crater ellipse angle w.r.t. x-direction in image plane
                batch_size : int
                    Return single detection feature, or create a batch for array-values
                convert_to_radians : bool
                    Whether to convert psi to radians inside method (default: True)
                max_iter : int
                    Maximum iterations (default: 10000)
                sort_ij : bool
                    Whether to sort triad features with I_ij being the lowest absolute value
        
                Yields
                ------
                crater_triad : np.ndarray
                    Triad indices (1x3)
                CoplanarInvariants
                     Associated CoplanarInvariants instance
        
                References
                ----------
                .. [1] Christian, J. A., Derksen, H., & Watkins, R. (2020). Lunar Crater Identification in Digital Images. https://arxiv.org/abs/2009.01228
                .. [2] Arnas, D., Fialho, M. A. A., & Mortari, D. (2017). Fast and robust kernel generators for star trackers. Acta Astronautica, 134 (August 2016), 291–302. https://doi.org/10.1016/j.actaastro.2017.02.016
                
        """
    def __init__(self, crater_triads, A_i, A_j, A_k, normalize_det = True):
        """
        Generates projective invariants [1] assuming craters are coplanar. Input is an array of crater matrices
                such as those generated using L{conic_matrix}.
        
                Parameters
                ----------
                crater_triads : np.ndarray
                    Crater triad indices (nx3) for slicing arrays
                A_i : np.ndarray
                    Crater representation first crater in triad
                A_j : np.ndarray
                    Crater representation second crater in triad
                A_k : np.ndarray
                    Crater representation third crater in triad
                normalize_det : bool
                    Set to True to normalize matrices to achieve det(A) = 1
        
                References
                ----------
                .. [1] Christian, J. A., Derksen, H., & Watkins, R. (2020). Lunar Crater Identification in Digital Images. https://arxiv.org/abs/2009.01228
                
        """
    def __len__(self):
        ...
    def get_pattern(self, permutation_invariant = False):
        """
        Get matching pattern using either permutation invariant features (eq. 134 from [1]) or raw projective
                invariants (p. 61 from [1]).
        
                Parameters
                ----------
                permutation_invariant : bool
                    Set this to True if permutation_invariants are needed.
        
                Returns
                -------
                np.ndarray
                    Array of features linked to the crater triads generated during initialisation.
        
                References
                ----------
                .. [1] Christian, J. A., Derksen, H., & Watkins, R. (2020). Lunar Crater Identification in Digital Images. https://arxiv.org/abs/2009.01228
                
        """
class PermutationInvariant:
    """
    
        Namespace for permutation invariants functions
        
    """
    @staticmethod
    def F1(x, y, z):
        ...
    @staticmethod
    def F2(x, y, z):
        ...
    @staticmethod
    def F3(x, y, z):
        ...
    @staticmethod
    def G1(x1, y1, z1, x2, y2, z2):
        ...
    @staticmethod
    def G2(x1, y1, z1, x2, y2, z2):
        ...
    @classmethod
    def F(cls, x, y, z):
        """
        Three-pair cyclic permutation invariant function.
        
                Parameters
                ----------
                x, y, z : int or float or np.ndarray
                    Values to generate cyclic permutation invariant features for
        
                Returns
                -------
                np.ndarray
                    Array containing cyclic permutation invariants F1, F2, F3
        
                References
                ----------
                .. [1] Christian, J. A., Derksen, H., & Watkins, R. (2020). Lunar Crater Identification in Digital Images. https://arxiv.org/abs/2009.01228
                
        """
    @classmethod
    def G(cls, x1, y1, z1, x2, y2, z2):
        ...
    @classmethod
    def G_tilde(cls, x1, y1, z1, x2, y2, z2):
        """
        
        
                Parameters
                ----------
                x1, y1, z1 : int or float or np.ndarray
                    First set of values to generate cyclic permutation invariant features for
                x2, y2, z2 : int or float or np.ndarray
                    Second set of values to generate cyclic permutation invariant features for
        
                Returns
                -------
                np.ndarray
                    Array containing cyclic permutation invariants G1, G2
        
                References
                ----------
                .. [1] Christian, J. A., Derksen, H., & Watkins, R. (2020). Lunar Crater Identification in Digital Images. https://arxiv.org/abs/2009.01228
        
                
        """
