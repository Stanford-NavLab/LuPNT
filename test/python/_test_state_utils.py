# def test_StateUtils():
#     # ELFO 9750.5 0.7 63.5 90 0 0
#     # LLO 100 0 28.5 0 0 0
#     # PCO 3000 0 75 0 90 0
#     a = 970
#     e = 0.6
#     i = np.deg2rad(65.5)
#     Omega = np.deg2rad(90.0)
#     w = np.deg2rad(0.0)
#     M = np.deg2rad(0.0)

#     coe_state = pnt.ClassicalOE([a, e, i, Omega, w, M], pnt.CoordSystem.MI)

#     mu = pnt.MU_MOON
#     coe_vec = coe_state.get_vector()

#     np.testing.assert_array_almost_equal(coe_vec, [a, e, i, Omega, w, M])

#     cart_vec_gmat = gmat.StateConversionUtil.KeplerianToCartesian(
#         mu, gmat.Rvector6(*coe_vec), "MA"
#     )
#     cart_vec_gmat = np.array([cart_vec_gmat.Get(i) for i in range(6)])
#     cart_vec_gmat[2:6] = np.deg2rad(cart_vec_gmat[2:6])

#     cart_vec = pnt.coe_to_cart(coe_vec, mu)
#     # np.testing.assert_array_almost_equal(cart_vec_gmat, cart_vec)
