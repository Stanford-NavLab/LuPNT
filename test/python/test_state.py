import pylupnt as pnt
import numpy as np
import pytest

try:
    from ..gmat.utils import data
except ImportError:
    from test.gmat.utils import data


class TestOrbitState:
    def test_OrbitStateRepres(self):
        attriibutes = [
            "CARTESIAN",
            "CLASSICAL_OE",
            "QUASI_NONSINGULAR_OE",
            "EQUINOTICAL_OE",
            "SINGULAR_ROE",
            "QUASINONSINGULAR_ROE",
            "DELAUNAY_OE",
        ]
        for attr in attriibutes:
            assert hasattr(pnt.OrbitStateRepres, attr)

    def test_ClassicalOE(self):
        # Constructor
        coe = data.coe_array_elfo
        coe_state = pnt.ClassicalOE(coe, pnt.CoordSystem.MI)
        attributes = ["a", "e", "i", "Omega", "w", "M"]

        # Getters
        for i, attr in enumerate(attributes):
            assert getattr(coe_state, attr) == coe[i]
        assert coe_state.coord_sys == pnt.CoordSystem.MI
        assert coe_state.state_repres == pnt.OrbitStateRepres.CLASSICAL_OE
        assert coe_state.size == 6
        assert np.allclose(coe_state.vector, coe)

        # Setters
        coe = data.coe_array_llo
        for i, attr in enumerate(attributes):
            setattr(coe_state, attr, coe[i])
            assert getattr(coe_state, attr) == coe[i]

        coe = data.coe_array_elfo
        coe_state.vector = coe
        assert np.allclose(coe_state.vector, coe)

        # Other
        assert (
            str(coe_state)
            == "<pylupnt.ClassicalOE "
            + np.array2string(coe, **data.array2string_kwargs)
            + ">"
        )

    def test_CartesianOrbitState(self):
        # Constructor
        rv = pnt.classical_to_cartesian(data.coe_array_elfo, pnt.MU_MOON)
        cart_state = pnt.CartesianOrbitState(rv, pnt.CoordSystem.MI)

        # Getters
        assert np.allclose(cart_state.r, rv[0:3])
        assert np.allclose(cart_state.v, rv[3:6])

        assert cart_state.coord_sys == pnt.CoordSystem.MI
        assert cart_state.state_repres == pnt.OrbitStateRepres.CARTESIAN
        assert cart_state.size == 6
        assert np.allclose(cart_state.vector, rv)

        # Setters
        rv = pnt.classical_to_cartesian(data.coe_array_llo, pnt.MU_MOON)
        cart_state.r = rv[0:3]
        cart_state.v = rv[3:6]
        assert np.allclose(cart_state.r, rv[0:3])
        assert np.allclose(cart_state.v, rv[3:6])

        rv = pnt.classical_to_cartesian(data.coe_array_elfo, pnt.MU_MOON)
        cart_state.vector = rv
        assert np.allclose(cart_state.vector, rv)

        # Other
        assert (
            str(cart_state)
            == "<pylupnt.CartesianOrbitState "
            + np.array2string(rv, **data.array2string_kwargs)
            + ">"
        )

    def test_QuasiNonsingularOE(self):
        # Constructor
        qns_oe = pnt.classical_to_quasi_nonsingular(data.coe_array_elfo)
        qns_oe_state = pnt.QuasiNonsingularOE(qns_oe, pnt.CoordSystem.MI)
        attributes = ["a", "u", "ex", "ey", "i", "Omega"]

        # Getters
        for i, attr in enumerate(attributes):
            assert getattr(qns_oe_state, attr) == qns_oe[i]
        assert qns_oe_state.coord_sys == pnt.CoordSystem.MI
        assert qns_oe_state.state_repres == pnt.OrbitStateRepres.QUASI_NONSINGULAR_OE
        assert qns_oe_state.size == 6
        assert np.allclose(qns_oe_state.vector, qns_oe)

        # Setters
        qns_oe = pnt.classical_to_quasi_nonsingular(data.coe_array_llo)
        for i, attr in enumerate(attributes):
            setattr(qns_oe_state, attr, qns_oe[i])
            assert getattr(qns_oe_state, attr) == qns_oe[i]

        qns_oe = pnt.classical_to_quasi_nonsingular(data.coe_array_elfo)
        qns_oe_state.vector = qns_oe
        assert np.allclose(qns_oe_state.vector, qns_oe)

        # Other
        assert (
            str(qns_oe_state)
            == "<pylupnt.QuasiNonsingularOE "
            + np.array2string(qns_oe, **data.array2string_kwargs)
            + ">"
        )

    def test_EquinoctialOE(self):
        # Constructor
        eq_oe = pnt.classical_to_equinoctial(data.coe_array_elfo)
        eq_oe_state = pnt.EquinoctialOE(eq_oe, pnt.CoordSystem.MI)
        attributes = ["a", "h", "k", "p", "q", "lon"]

        # Getters
        for i, attr in enumerate(attributes):
            assert getattr(eq_oe_state, attr) == eq_oe[i]
        assert eq_oe_state.coord_sys == pnt.CoordSystem.MI
        assert eq_oe_state.state_repres == pnt.OrbitStateRepres.EQUINOTICAL_OE
        assert eq_oe_state.size == 6

        # Setters
        eq_oe = pnt.classical_to_equinoctial(data.coe_array_llo)
        for i, attr in enumerate(attributes):
            setattr(eq_oe_state, attr, eq_oe[i])
            assert getattr(eq_oe_state, attr) == eq_oe[i]

        eq_oe = pnt.classical_to_equinoctial(data.coe_array_elfo)
        eq_oe_state.vector = eq_oe
        assert np.allclose(eq_oe_state.vector, eq_oe)

        # Other
        assert (
            str(eq_oe_state)
            == "<pylupnt.EquinoctialOE "
            + np.array2string(eq_oe, **data.array2string_kwargs)
            + ">"
        )

    def test_SingularROE(self):
        # Constructor
        roe = data.roe_array_1
        roe_state = pnt.SingularROE(roe, pnt.CoordSystem.MI)
        attributes = ["ada", "adM", "ade", "adw", "adi", "adOmega"]

        # Getters
        for i, attr in enumerate(attributes):
            assert getattr(roe_state, attr) == roe[i]
        assert roe_state.coord_sys == pnt.CoordSystem.MI
        assert roe_state.state_repres == pnt.OrbitStateRepres.SINGULAR_ROE
        assert roe_state.size == 6

        # Setters
        roe = data.roe_array_2
        for i, attr in enumerate(attributes):
            setattr(roe_state, attr, roe[i])
            assert getattr(roe_state, attr) == roe[i]

        roe = data.roe_array_1
        roe_state.vector = roe
        assert np.allclose(roe_state.vector, roe)

        # Other
        assert (
            str(roe_state)
            == "<pylupnt.SingularROE "
            + np.array2string(roe, **data.array2string_kwargs)
            + ">"
        )

    def test_QuasiNonsingularROE(self):
        # Constructor
        roe = data.roe_array_1
        roe_state = pnt.QuasiNonsingularROE(roe, pnt.CoordSystem.MI)
        attributes = ["ada", "adl", "adex", "adey", "adix", "adiy"]

        # Getters
        for i, attr in enumerate(attributes):
            assert getattr(roe_state, attr) == roe[i]
        assert roe_state.coord_sys == pnt.CoordSystem.MI
        assert roe_state.state_repres == pnt.OrbitStateRepres.QUASINONSINGULAR_ROE
        assert roe_state.size == 6

        # Setters
        roe = data.roe_array_2
        for i, attr in enumerate(attributes):
            setattr(roe_state, attr, roe[i])
            assert getattr(roe_state, attr) == roe[i]

        roe = data.roe_array_1
        roe_state.vector = roe
        assert np.allclose(roe_state.vector, roe)

        # Other
        assert (
            str(roe_state)
            == "<pylupnt.QuasiNonsingularROE "
            + np.array2string(roe, **data.array2string_kwargs)
            + ">"
        )


if __name__ == "__main__":
    pytest.main([__file__])
