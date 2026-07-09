Estimation Filters
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT's estimation
filters: the recursive Kalman family (``EKF``, ``UKF``, ``UDUEKF``, ``SRIF``,
``SchmidtEKF``), the batch weighted least-squares estimator, and the
process-noise models that drive them.  The stochastic-cloning UD filter and
smoother follow the author's PhD thesis, Chapter 6.2; the corresponding
equations are cited inline.  The concrete classes live under
``cpp/lupnt/numerics/filters/`` and ``cpp/lupnt/dynamics/clock_dynamics.*``.

State, Covariance, and Callback Contract
-------------------------------------------------------------------

Every recursive filter derives from ``Filter`` and maintains a state estimate
and covariance

.. math::

   \hat{x} \in \mathbb{R}^{n},
   \qquad
   P \in \mathbb{R}^{n\times n},
   \qquad P = P^\mathsf{T} \succeq 0 ,

with prior (pre-update) and posterior (post-update) copies
:math:`(\hat{x}^-, P^-)` and :math:`(\hat{x}^+, P^+)`.  The filter is driven by
three registered callbacks:

.. math::

   \texttt{f\_dyn\_}: (x, t_0, t_f, u) \mapsto \big(x(t_f),\ F\big),
   \qquad
   F = \frac{\partial x(t_f)}{\partial x(t_0)} ,

.. math::

   \texttt{f\_proc\_}: (x, t_0, t_f) \mapsto Q,
   \qquad
   \texttt{f\_meas\_}: x \mapsto \big(\hat{z},\ H,\ R\big),
   \quad H = \frac{\partial z}{\partial x} .

The estimation loop alternates ``Predict(t)`` (propagate to :math:`t`) and
``Update(z)`` (assimilate a measurement vector).  The state-transition matrix
:math:`F` and measurement Jacobian :math:`H` are supplied by
``f_dyn_``/``f_meas_`` (produced from plain models via
``GetFilterDynamicsFunction`` / ``GetFilterMeasurementFunction`` using
autodiff).

The callback typedefs and the shared ``x_`` / ``P_`` / ``F_`` / ``H_`` / ``R_``
members are declared on the base class ``cpp/lupnt/numerics/filters/filter.h``.

Extended Kalman Filter (EKF)
-------------------------------------------------------------------

**Predict.**  ``EKF::Predict`` evaluates the process noise at the prior state,
propagates the state and STM, and forms the propagated covariance:

.. math::

   Q = \texttt{f\_proc\_}(\hat{x}, t_0, t_f),
   \qquad
   \hat{x}^- = \texttt{f\_dyn\_}(\hat{x}, t_0, t_f, u,\, F),

.. math::

   \bar{P} = F P F^\mathsf{T},
   \qquad
   P^- = \bar{P} + Q .

When a process-noise mapping :math:`G` is registered
(``SetProcessNoiseMappingMatrix``) the noise is mapped first,
:math:`Q \leftarrow G Q G^\mathsf{T}`.  The propagated-only block
:math:`\bar{P} = F P F^\mathsf{T}` is retained (``GetCovarianceBar``) for
adaptive process-noise estimation.

Implemented by ``cpp/lupnt/numerics/filters/ekf.cc :: EKF::Predict``:

.. code-block:: cpp

   Q_ = f_proc_(x_, t_, t);
   if (use_process_noise_mapping_) Q_ = G_ * Q_ * G_.transpose();
   x_ = f_dyn_(x_, t_, t, u, &F_);
   P_bar_ = F_ * P_ * F_.transpose();
   P_ = P_bar_ + Q_;

**Update.**  ``EKF::Update`` linearizes the measurement, forms the innovation
covariance and residual, rejects outliers, and applies a Joseph-form update:

.. math::

   \hat{z} = \texttt{f\_meas\_}(\hat{x}^-,\, H,\, R),
   \qquad
   S = R + H P^- H^\mathsf{T},
   \qquad
   \delta z = z - \hat{z} ,

.. math::

   K = P^- H^\mathsf{T} S^{-1},
   \qquad
   \delta x = K\,\delta z,
   \qquad
   \hat{x}^+ = \hat{x}^- + \delta x ,

.. math::

   P^+ = (I - K H)\,P^-\,(I - K H)^\mathsf{T} + K R K^\mathsf{T} .

The state-correction covariance :math:`\Sigma_{\delta x} = K S K^\mathsf{T}`
is recorded.

Implemented by ``cpp/lupnt/numerics/filters/ekf.cc :: EKF::Update``:

.. code-block:: cpp

   z_prior_ = f_meas_(x_, &H_, &R_);
   S_ = R_ + H_ * P_ * H_.transpose();
   dz_ = z_true_ - z_prior_;
   K_ = P_ * H_.transpose() * S_.inverse();
   dx_ = K_ * dz_;  x_ = x_ + dx_;
   MatXd G = I - K_ * H_;
   P_ = G * P_ * G.transpose() + K_ * R_ * K_.transpose();  // Joseph form

**Outlier rejection.**  Before the gain is applied, ``RemoveOutliers`` drops
any measurement component whose normalized residual exceeds the threshold
:math:`\tau` (default :math:`3\sigma`),

.. math::

   \frac{|\delta z_i|}{\sqrt{S_{ii}}} > \tau
   \;\Rightarrow\;
   \text{remove } i ,

shrinking :math:`\delta z, H, R` and recomputing
:math:`S = H P^- H^\mathsf{T} + R`.  A custom
``FaultDetectionFunction`` may replace this default.

Implemented by ``cpp/lupnt/numerics/filters/ekf.cc :: EKF::RemoveOutliers``:

.. code-block:: cpp

   VecXd ratio = dz_.array().abs() / S_.diagonal().array().sqrt();
   VecXb is_outlier = ratio.array() > outlier_threshold_;
   // ... drop flagged rows from dz_, H_, R_ ...
   S_ = H_ * P_prior_ * H_.transpose() + R_;

**Schmidt / consider states.**  ``SetConsiderStateCount(n_c)`` (also exposed
as ``SchmidtEKF``) marks the trailing :math:`n_c` states as consider
parameters.  They enter :math:`F`, :math:`H`, :math:`K`, and the covariance
bookkeeping, but their mean is never corrected: the last :math:`n_c` rows of
:math:`K` are zeroed before :math:`\delta x = K\,\delta z`.  The Joseph-form
update remains valid for this non-optimal gain, so their (cross-)covariance
still updates consistently.

Implemented by the ``n_consider_`` branch of
``cpp/lupnt/numerics/filters/ekf.cc :: EKF::Update``
(``SchmidtEKF`` is a thin subclass in ``schmidt_ekf.h``):

.. code-block:: cpp

   if (n_consider_ > 0) K_.bottomRows(n_consider_).setZero();  // never correct consider states
   dx_ = K_ * dz_;

**RTS smoother.**  With the logged prior/posterior states, covariances, and
STMs, ``EKF::UpdateSmoother`` runs the standard Rauch-Tung-Striebel backward
recursion:

.. math::

   A_k = P_{k|k}\,F_{k+1}^\mathsf{T}\,P_{k+1|k}^{-1} ,

.. math::

   \hat{x}_{k|N} = \hat{x}_{k|k} + A_k\big(\hat{x}_{k+1|N} - \hat{x}_{k+1|k}\big),

.. math::

   P_{k|N} = P_{k|k} + A_k\big(P_{k+1|N} - P_{k+1|k}\big)A_k^\mathsf{T} .

Implemented by ``cpp/lupnt/numerics/filters/ekf.cc :: EKF::UpdateSmoother``:

.. code-block:: cpp

   MatXd A_k = P_kk * stm_k.transpose() * P_k1k.inverse();
   x_sm_[k] = x_kk + A_k * (x_sm_[k + 1] - x_k1k);
   P_sm_[k] = P_kk + A_k * (P_sm_[k + 1] - P_k1k) * A_k.transpose();

Unscented Kalman Filter (UKF)
-------------------------------------------------------------------

``UKF`` replaces linearization with the unscented transform.  For an
:math:`n`-dimensional state it uses :math:`2n+1` sigma points with scaling

.. math::

   \lambda = \alpha^2 (n + \kappa) - n ,

(defaults :math:`\alpha = 10^{-3}`, :math:`\beta = 2`, :math:`\kappa = 0`) and
weights

.. math::

   w^m_0 = \frac{\lambda}{n+\lambda},
   \qquad
   w^c_0 = \frac{\lambda}{n+\lambda} + (1 - \alpha^2 + \beta),
   \qquad
   w^m_i = w^c_i = \frac{1}{2(n+\lambda)} .

The sigma points about :math:`(\hat{x}, P)` are

.. math::

   \mathcal{X}_0 = \hat{x},
   \qquad
   \mathcal{X}_i = \hat{x} \pm \Big[\sqrt{(n+\lambda)\,P}\Big]_i ,

with the matrix square root from a Cholesky factor.  **Predict** propagates
each sigma point through ``f_dyn_`` and recombines with process noise:

.. math::

   \hat{x}^- = \sum_i w^m_i\,\mathcal{X}_i^{-},
   \qquad
   P^- = \sum_i w^c_i (\mathcal{X}_i^{-} - \hat{x}^-)(\mathcal{X}_i^{-} - \hat{x}^-)^\mathsf{T} + Q .

**Update** transforms the sigma points through ``f_meas_`` to measurement
points :math:`\mathcal{Z}_i`, and forms

.. math::

   \hat{z} = \sum_i w^m_i \mathcal{Z}_i,
   \qquad
   S = \sum_i w^c_i (\mathcal{Z}_i - \hat{z})(\mathcal{Z}_i - \hat{z})^\mathsf{T} + R ,

.. math::

   P_{xz} = \sum_i w^c_i (\mathcal{X}_i - \hat{x}^-)(\mathcal{Z}_i - \hat{z})^\mathsf{T},
   \qquad
   K = P_{xz} S^{-1} ,

.. math::

   \hat{x}^+ = \hat{x}^- + K(z - \hat{z}),
   \qquad
   P^+ = P^- - K S K^\mathsf{T} .

The weights are set in
``cpp/lupnt/numerics/filters/ukf.cc :: UKF::InitializeUkfParams``; the sigma
points in ``UKF::ComputeSigmaPoints``; predict/update in ``UKF::Predict`` /
``UKF::Update`` (all in ``ukf.cc``, declared in ``ukf.h``):

.. code-block:: cpp

   // InitializeUkfParams
   lambda_  = alpha_ * alpha_ * (n_x_ + kappa_) - n_x_;
   w_m_(0)  = lambda_ / (n_x_ + lambda_);
   w_c_(0)  = lambda_ / (n_x_ + lambda_) + (1.0 - alpha_ * alpha_ + beta_);

   // ComputeSigmaPoints -- columns of chol((n_x_ + lambda_) * cov)
   MatXd chol_cov = ((n_x_ + lambda_) * cov).llt().matrixL();
   sigma_points.col(i + 1)        = state + chol_cov.col(i);
   sigma_points.col(i + 1 + n_x_) = state - chol_cov.col(i);

UD-Factorized Filter (Bierman-Thornton)
-------------------------------------------------------------------

``UDUEKF`` maintains the covariance in factored form

.. math::

   P = U\,\operatorname{diag}(D)\,U^\mathsf{T} ,

with :math:`U` unit upper-triangular and :math:`D` diagonal (positive), which
guarantees symmetry and positive semi-definiteness by construction.  This
mirrors the thesis UD factorization :math:`\hat{P}_{k|m} = U_{k|m} D_{k|m}
U_{k|m}^\mathsf{T}` (Iiyama, Ch. 6.2.1, Eq. 6.17; Algorithm 1).

**Time update (Thornton MWGS).**  The predict step never re-forms :math:`P`.
With the augmented coefficient and weight matrices

.. math::

   Y = \big[\,F U \;\big|\; G\,\big],
   \qquad
   \tilde{D} = \operatorname{blkdiag}(D,\ Q) ,

a modified weighted Gram-Schmidt orthogonalization returns
:math:`(D^-, U^-)` such that

.. math::

   U^- \operatorname{diag}(D^-) (U^-)^\mathsf{T}
   = Y\,\tilde{D}\,Y^\mathsf{T}
   = F P F^\mathsf{T} + G Q G^\mathsf{T} ,

computed without forming the product (thesis Ch. 6.2.2, Eqs. 6.20-6.23).  The
process noise :math:`Q` must be diagonal; correlated noise is supplied through
the mapping matrix :math:`G`.

Implemented by ``cpp/lupnt/numerics/filters/udu_filter.cc :: UDUEKF::Predict``
(the MWGS itself is ``ModifiedGramSchmidt`` in ``udu_utils.cc``):

.. code-block:: cpp

   MatXd Y = MatXd::Zero(n, n + m);
   Y.leftCols(n)  = F_ * U_;
   Y.rightCols(m) = G_;
   MatXd D_tilde = MatXd::Zero(n + m, n + m);
   D_tilde.topLeftCorner(n, n)     = D_diag_.asDiagonal();
   D_tilde.bottomRightCorner(m, m) = Q_;
   std::tie(D_diag_, U_) = ModifiedGramSchmidt(D_tilde, Y);

**Measurement update (Carlson rank-one).**  Measurements are processed
sequentially, one scalar component at a time (requiring diagonal :math:`R`).
For a scalar row :math:`h` with variance :math:`r`, define

.. math::

   f = U^\mathsf{T} h,
   \qquad
   v = D\,f,
   \qquad
   \alpha_0 = r,
   \qquad
   \alpha_k = \alpha_{k-1} + v_k f_k ,

.. math::

   D^+_k = \frac{\alpha_{k-1}}{\alpha_k}\,D_k ,

with the Carlson recursion updating the columns of :math:`U` in place and
accumulating the gain :math:`K = \beta / \alpha_n` (``UDUEKF::CarlsonUpdate``).
This is the rank-one UD downdate of the thesis (Ch. 6.2.3, Eqs. 6.24-6.27;
Algorithm 2).

Implemented by
``cpp/lupnt/numerics/filters/udu_filter.cc :: UDUEKF::CarlsonUpdate``:

.. code-block:: cpp

   VecXd f = U_.transpose() * h;         // f = U^T h
   VecXd v = D_diag_.asDiagonal() * f;   // v = D f
   alpha(0) = r;
   for (int k = 0; k < n; ++k) {
     alpha(k + 1) = alpha(k) + v(k) * f(k);
     D_plus(k)    = alpha(k) / alpha(k + 1) * D_diag_(k);
     // ... rank-one column update of U ...
   }
   VecXd k_gain = beta / alpha(n);

Square-Root Information Filter (SRIF)
-------------------------------------------------------------------

``SRIF`` (a subclass of ``EKF``) stores the covariance as an
upper-triangular information square root :math:`R_\text{info}`:

.. math::

   R_\text{info}^\mathsf{T} R_\text{info} = P^{-1},
   \qquad
   P = R_\text{info}^{-1} R_\text{info}^{-\mathsf{T}} .

The dense :math:`P`, :math:`\hat{x}` are kept in sync after each step so the
inherited outlier rejection, fault detection, and RTS smoother work unchanged.

**Time update (Dyer-McReynolds).**  With the process-noise square root
:math:`S` (:math:`Q = S S^\mathsf{T}`) and
:math:`R_\Phi = R_\text{info}\,\Phi^{-1}`, a Householder QR triangularization
of the array

.. math::

   \begin{bmatrix}
     I & 0 \\
     -R_\Phi S & R_\Phi
   \end{bmatrix}
   \;\xrightarrow{\ \text{QR}\ }\;
   \begin{bmatrix}
     * & * \\
     0 & \bar{R}_\text{info}
   \end{bmatrix}

zeros the bottom-left block; the bottom-right block is the a-priori
information square root at :math:`t_f`.  The estimate is re-centered on
:math:`\hat{x}` each step, so the data column is zero.

Implemented by ``cpp/lupnt/numerics/filters/srif.cc :: SRIF::Predict``
(``QrRFactor`` is the file-local Householder helper):

.. code-block:: cpp

   MatXd RPhiInv = R_info_ * F_.inverse();
   MatXd A = MatXd::Zero(2 * n, 2 * n);
   A.topLeftCorner(n, n)     = MatXd::Identity(n, n);
   A.bottomLeftCorner(n, n)  = -RPhiInv * S;   // Q = S S^T
   A.bottomRightCorner(n, n) = RPhiInv;
   R_info_ = QrRFactor(A).bottomRightCorner(n, n);

**Measurement update.**  The rows are whitened by :math:`R^{-1/2}` (with
:math:`R = L L^\mathsf{T}`, :math:`H_w = L^{-1}H`,
:math:`\delta z_w = L^{-1}\delta z`) and folded into the information array by
a QR triangularization

.. math::

   \begin{bmatrix}
     R_\text{info} & 0 \\
     H_w & \delta z_w
   \end{bmatrix}
   \;\xrightarrow{\ \text{QR}\ }\;
   \begin{bmatrix}
     \hat{R}_\text{info} & \hat{b} \\
     0 & *
   \end{bmatrix},
   \qquad
   \delta x = \hat{R}_\text{info}^{-1}\hat{b},
   \qquad
   \hat{x}^+ = \hat{x}^- + \delta x .

Implemented by
``cpp/lupnt/numerics/filters/srif.cc :: SRIF::InformationUpdate``:

.. code-block:: cpp

   MatXd L = llt.matrixL();                      // R = L L^T
   MatXd H_w  = L.triangularView<Eigen::Lower>().solve(H_);
   VecXd dz_w = L.triangularView<Eigen::Lower>().solve(dz_);
   // QR of [R_info_ 0; H_w dz_w] -> R_hat (top-left), b_hat (top of last col)
   dx_ = R_info_.triangularView<Eigen::Upper>().solve(b_hat);
   x_  = x_ + dx_;

Batch Weighted Least Squares
-------------------------------------------------------------------

``RunBatchFilter`` solves the nonlinear weighted least-squares problem by
Gauss-Newton iteration.  At each iteration it stacks all measurement
residuals :math:`\delta z = z - h(\hat{x})`, Jacobians :math:`H`, and weights
:math:`W = \operatorname{diag}(1/\sigma_i^2)`, and solves the normal equations

.. math::

   \big(H^\mathsf{T} W H\big)\,\delta x = H^\mathsf{T} W\,\delta z,
   \qquad
   \hat{x} \leftarrow \hat{x} + \delta x ,

via an LDLT solve (``SolveWeightedLeastSquares``).  When initialization
constraints are enabled, a soft prior is appended by augmenting

.. math::

   H \leftarrow \begin{bmatrix} H \\ I \end{bmatrix},
   \qquad
   \delta z \leftarrow \begin{bmatrix} \delta z \\ x_0 - \hat{x} \end{bmatrix},
   \qquad
   W \leftarrow \operatorname{blkdiag}\!\big(W,\ \operatorname{diag}(1/\sigma_{0}^2)\big) .

Iteration stops when :math:`\|\delta x\| <` ``convergence_tol`` or the
iteration limit is reached.  The final covariance is the inverse information
(measurement Jacobians only, via SVD pseudo-inverse for rank-deficient cases)

.. math::

   P = \big(H^\mathsf{T} W H\big)^{-1},
   \qquad
   \text{PDOP} = \sqrt{\operatorname{tr}\!\big(P_{[0:3,0:3]}\big)} .

Implemented by
``cpp/lupnt/numerics/filters/batch_filter.cc :: RunBatchFilter``, with the
normal-equation solve in ``SolveWeightedLeastSquares``, the soft prior in
``AddInitializationConstraints``, and the PDOP in ``CalculatePDOP`` (all in
``batch_filter.cc``):

.. code-block:: cpp

   // SolveWeightedLeastSquares
   MatXd HTW  = H.transpose() * weights.asDiagonal();
   return (HTW * H).ldlt().solve(HTW * residuals);   // (H^T W H) dx = H^T W r

   // final covariance -- SVD pseudo-inverse of the measurement information
   Eigen::JacobiSVD<MatXd> svd(information_matrix, ComputeFullU | ComputeFullV);
   covariance_matrix = svd.solve(MatXd::Identity(n_state, n_state));

Process-Noise Models
-------------------------------------------------------------------

**Continuous white-noise acceleration (CWNA).**  For a Cartesian
:math:`[r\ (3);\ v\ (3)]` state over a step :math:`\Delta t`, with per-axis
acceleration PSD :math:`q_a` [:math:`\mathrm{m}^2/\mathrm{s}^3`],
``CwnaProcessNoise(dt, accel_psd)`` returns the :math:`6\times 6` block

.. math::

   Q = q_a
   \begin{bmatrix}
     \tfrac{\Delta t^3}{3}\, I_3 & \tfrac{\Delta t^2}{2}\, I_3 \\[4pt]
     \tfrac{\Delta t^2}{2}\, I_3 & \Delta t\, I_3
   \end{bmatrix}.

Implemented by
``cpp/lupnt/numerics/filters/srif.cc :: CwnaProcessNoise``:

.. code-block:: cpp

   const double q3 = accel_psd * dt * dt * dt / 3.0;
   const double q2 = accel_psd * dt * dt / 2.0;
   const double q1 = accel_psd * dt;
   for (int k = 0; k < 3; ++k) {
     Q(k, k) = q3;  Q(k, k + 3) = q2;  Q(k + 3, k) = q2;  Q(k + 3, k + 3) = q1;
   }

**Clock two- and three-state noise.**  ``ClockDynamics`` propagates a clock
error state with the polynomial transition matrices

.. math::

   \Phi_2(\Delta t) = \begin{bmatrix} 1 & \Delta t \\ 0 & 1 \end{bmatrix},
   \qquad
   \Phi_3(\Delta t) = \begin{bmatrix}
     1 & \Delta t & \tfrac{\Delta t^2}{2} \\
     0 & 1 & \Delta t \\
     0 & 0 & 1
   \end{bmatrix},

and closed-form random-walk process-noise covariances built from the
oscillator PSD coefficients :math:`(q_1, q_2, q_3)` (white-frequency,
random-walk-frequency, and frequency-drift; ``GetClockValues``).  For the
two-state model (``TwoStateNoise``):

.. math::

   Q_2 =
   \begin{bmatrix}
     q_1 \Delta t + q_2\dfrac{\Delta t^3}{3} & q_2\dfrac{\Delta t^2}{2} \\[8pt]
     q_2\dfrac{\Delta t^2}{2} & q_2 \Delta t
   \end{bmatrix},

and for the three-state model (``ThreeStateNoise``):

.. math::

   Q_3 =
   \begin{bmatrix}
     q_1\Delta t + q_2\dfrac{\Delta t^3}{3} + q_3\dfrac{\Delta t^5}{20}
       & q_2\dfrac{\Delta t^2}{2} + q_3\dfrac{\Delta t^4}{8}
       & q_3\dfrac{\Delta t^3}{6} \\[8pt]
     q_2\dfrac{\Delta t^2}{2} + q_3\dfrac{\Delta t^4}{8}
       & q_2\Delta t + q_3\dfrac{\Delta t^3}{3}
       & q_3\dfrac{\Delta t^2}{2} \\[8pt]
     q_3\dfrac{\Delta t^3}{6}
       & q_3\dfrac{\Delta t^2}{2}
       & q_3\Delta t
   \end{bmatrix}.

These are scaled to the configured clock-bias unit (seconds, meters,
kilometers) by :math:`Q \leftarrow \sigma^2 Q` with
:math:`\sigma = \texttt{SecondsToBiasUnitScale}` (:math:`1`, :math:`c`, or
:math:`c\cdot\mathrm{km/m}`).  This matches the coupled orbit-clock process
noise of the thesis (Ch. 6.1.2, Eq. 6.6), whose diffusion coefficients
:math:`q_1, q_2, q_3` correspond to phase, frequency, and frequency-drift
(aging) noise.

Implemented in ``cpp/lupnt/dynamics/clock_dynamics.cc``:
``ClockDynamics::GetClockValues`` (the per-oscillator PSDs),
``ClockDynamics::TwoStatePhi`` / ``ThreeStatePhi`` (the :math:`\Phi` above), and
``ClockDynamics::TwoStateNoise`` / ``ThreeStateNoise`` (the :math:`Q` above,
with the bias-unit overload applying :math:`\sigma^2`):

.. code-block:: cpp

   // ThreeStateNoise
   Q(0, 0) = q1 * dt + q2 * dt3 / 3 + q3 * dt5 / 20.0;
   Q(0, 1) = q2 * dt2 / 2.0 + q3 * dt4 / 8.0;   Q(0, 2) = q3 * dt3 / 6.0;
   Q(1, 1) = q2 * dt + q3 * dt3 / 3.0;          Q(1, 2) = q3 * dt2 / 2.0;
   Q(2, 2) = q3 * dt;   // (lower triangle mirrored)

   // bias-unit overload: Q_unit = sigma^2 Q,  sigma = SecondsToBiasUnitScale(unit)
   return ThreeStateNoise(clk_model, dt) * scale * scale;

**Adaptive process noise.**  ``ProcessNoise`` implements a sliding-window
adaptive scheme (``ProcessNoiseAlgorithm`` :math:`\in`
:math:`\{\text{SNC}, \text{ASNC}, \text{ADMC}\}`) that re-estimates the
acceleration-noise level from the recent history of state corrections
:math:`\delta x`, correction covariances :math:`\Sigma_{\delta x}`, and the
propagated/posterior covariances :math:`\bar{P}, P^+`, clamped to configured
:math:`[\sigma_\text{min}^2, \sigma_\text{max}^2]` limits.  It is consumed as
a ``ProcessNoiseFunction`` via ``GetProcessNoiseFunction``.

Implemented by ``cpp/lupnt/numerics/filters/adaptive_process_noise.cc`` (class
``ProcessNoise``, declared in ``adaptive_process_noise.h``); the per-step
history :math:`(\delta x, \Sigma_{\delta x}, \bar P, P^+)` is fed in through
``ProcessNoise::Update`` and the estimate is exposed via
``ProcessNoise::GetProcessNoiseFunction``.

Stochastic-Cloning UD Filter and Smoother
-------------------------------------------------------------------

For time-differenced carrier-phase (TDCP) measurements, which depend on both
the current and a previous state, LuPNT augments the state with a cloned copy
of the previous state (``UDUStochasticCloningEKF``) and maintains the full UD
factorization of the augmented covariance.  This follows the author's thesis,
Chapter 6.2 (Iiyama, *Lunar ODTS using GNSS*), whose contract is reproduced
here for reference.

The augmented state and its UD-factored, time-propagated covariance are
(thesis Eqs. 6.18-6.23)

.. math::

   \mathbf{x}^{\mathrm{sc}}_{k-1}
   = \begin{bmatrix} \hat{\mathbf{x}}_{k-1|k-1} \\ \hat{\mathbf{x}}_{k-1|k-1} \end{bmatrix},
   \qquad
   \hat{\mathbf{U}}^{\mathrm{sc}}_{k|k-1}
   = \begin{bmatrix}
       \mathbf{G}_k & \boldsymbol{\Phi}_{k,k-1}\hat{\mathbf{U}}_{k-1|k-1} \\
       \mathbf{0} & \hat{\mathbf{U}}_{k-1|k-1}
     \end{bmatrix},
   \qquad
   \hat{\mathbf{D}}^{\mathrm{sc}}_{k|k-1}
   = \begin{bmatrix}
       \mathbf{Q}^D_k & \mathbf{0} \\
       \mathbf{0} & \hat{\mathbf{D}}_{k-1|k-1}
     \end{bmatrix},

so that the top block is the propagated state (via the variational STM
:math:`\boldsymbol{\Phi}_{k,k-1}`) and the bottom block is the static clone.
In LuPNT the augmented STM is
:math:`F = \left[\begin{smallmatrix} F_\text{base} & 0 \\ I & 0 \end{smallmatrix}\right]`
and the Thornton MWGS time update acts on the current block only (process
noise mapped through the base :math:`G`).

Implemented by
``cpp/lupnt/numerics/filters/udu_filter.cc :: UDUStochasticCloningEKF::Predict``
(augmented STM assembly and the current-block-only MWGS):

.. code-block:: cpp

   F_.topLeftCorner(n, n)    = F_base;                 // [ F_base 0 ]
   F_.bottomLeftCorner(n, n) = MatXd::Identity(n, n);  // [ I      0 ]
   MatXd E = MatXd::Zero(n2, m);  E.topRows(n) = G_;   // noise on current block only
   MatXd Y = MatXd::Zero(n2, n2 + m);
   Y.leftCols(n2) = F_ * U_;  Y.rightCols(m) = E;
   std::tie(D_diag_, U_) = ModifiedGramSchmidt(D_tilde, Y);

The measurement Jacobian with respect to the augmented state splits into
current- and cloned-state parts (thesis Eq. 6.25),

.. math::

   \mathbf{H}_{k,m}
   = \begin{bmatrix} \mathbf{H}^{(1)}_{k,m} & \mathbf{H}^{(2)}_{k,m} \end{bmatrix},
   \qquad
   \mathbf{H}^{(1)}_{k,m} = \frac{\partial h}{\partial \mathbf{X}_k},
   \qquad
   \mathbf{H}^{(2)}_{k,m} = \frac{\partial h}{\partial \mathbf{X}_{k-1}} ,

with :math:`\mathbf{H}^{(2)} = \mathbf{0}` for current-state-only
measurements (e.g. pseudorange).  Measurements are assimilated sequentially by
the Carlson rank-one UD update (Eqs. 6.24-6.30).  Ordering delayed-state
(TDCP) measurements before current-state-only measurements lets the
current-state update skip the cloned UD sub-blocks
:math:`\hat{U}_{22}`, reducing cost (thesis Ch. 6.2.4).

**Fixed-interval smoothing.**  For post-processing, the standard RTS smoother
(thesis Eqs. 6.33-6.35)

.. math::

   \hat{\mathbf{x}}_{k|N} = \hat{\mathbf{x}}_{k|k}
     + \mathbf{A}_k\big(\hat{\mathbf{x}}_{k+1|N} - \hat{\mathbf{x}}_{k+1|k}\big),
   \qquad
   \mathbf{A}_k = \hat{\mathbf{P}}_{k|k}\boldsymbol{\Phi}_{k+1,k}^\mathsf{T}
     \hat{\mathbf{P}}_{k+1|k}^{-1} ,

relies on the Markov assumption (Eq. 6.36), which TDCP measurements violate
because :math:`\mathbf{z}_{k+1} = h_1(\mathbf{x}_{k+1}) - h_2(\mathbf{x}_k) +
\mathbf{v}_{k+1}` couples consecutive epochs (Eq. 6.37).  The thesis therefore
derives the smoother gain from the joint stochastic-cloning posterior
(Eqs. 6.38-6.41), giving

.. math::

   \hat{\mathbf{x}}_{k|N} = \hat{\mathbf{x}}_{k+1|k+1}
     + \mathbf{J}_k\big(\hat{\mathbf{x}}_{k+1|N} - \hat{\mathbf{x}}_{k+1|k+1}\big),
   \qquad
   \mathbf{J}_k = \hat{\mathbf{P}}_{k+1,k|k+1}^\mathsf{T}\,\hat{\mathbf{P}}_{k+1|k+1}^{-1} ,

where :math:`\hat{\mathbf{P}}_{k+1,k|k+1}` is the posterior cross-covariance
between the current and cloned states, which captures the
measurement-induced correlation that the standard RTS gain
:math:`\mathbf{A}_k` misses.  The gain is solved from the stored UD factors
without explicit inversion (Eq. 6.42) by backward substitution, diagonal
scaling, and forward substitution.

.. note::

   ``EKF::UpdateSmoother`` implements the **standard** RTS recursion
   (:math:`\mathbf{A}_k` above), which is exact when measurements depend only
   on the current state.  The cross-covariance-based smoother gain
   :math:`\mathbf{J}_k` is the thesis extension required when delayed-state
   (TDCP) measurements are present; it operates on the augmented
   stochastic-cloning covariance maintained by ``UDUStochasticCloningEKF``.

Model Boundaries
-------------------------------------------------------------------

* ``UDUEKF`` and ``UDUStochasticCloningEKF`` require diagonal process noise
  :math:`Q` and diagonal measurement noise :math:`R`; correlated process noise
  must be supplied through the mapping matrix :math:`G`.
* ``EKF``/``SRIF`` re-linearize about the current estimate each step (extended
  filters); the STM :math:`F` and Jacobian :math:`H` come from the registered
  callbacks.
* Outlier rejection and fault detection are shared across ``EKF``, ``SRIF``,
  and ``UDUEKF`` (SRIF/UDU reuse ``EKF::RemoveOutliers``).
* The batch estimator's reported covariance uses the measurement information
  only (initialization prior excluded).
