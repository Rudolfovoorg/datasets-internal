import time
import cvxpy as cvx
import numpy as np
import operator
import json
import sys


def qubo2maxcut(qubo: np.ndarray) -> np.ndarray:
    """Convert qubo to adjacency matrix that biqbin can read

    Args:
        qubo (np.ndarray): qubo as 2d numpy array
    Returns:
        np.ndarray: adjacency matrix for max cut problem
    """
    q_sym = 1/2*(qubo.T + qubo)
            
    Qe_plus_c = -np.array([(np.sum(q_sym, 1))])
    np.fill_diagonal(q_sym, 0)

    return np.block([
        [q_sym, Qe_plus_c.T],
        [Qe_plus_c, np.zeros((1, 1))]
    ])


def maxcut_sdp(MC, solver=cvx.MOSEK, verbose=False, **solver_kwargs):
    """
    Compute sdp relaxation of QUBO via MAX-CUT.
    solvers = [cvx.SCS, cvx.CVXOPT, cvx.MOSEK]
    MOSEK licence: https://www.mosek.com/products/version-11/
    """
    #MC = qubo2maxcut(Q) # Convert to a MAX-CUT 
    L = np.diag(MC.sum(axis=0))-MC # Compute a Laplacian matrix of G
    n, _ = L.shape
    
    # SDP solution
    X = cvx.Variable((n, n), PSD=True)
    obj = 0.25 * cvx.trace(L @ X)
    constr = [cvx.diag(X) == 1]
    problem = cvx.Problem(cvx.Maximize(obj), constraints=constr)
    problem.solve(solver=solver, verbose=verbose, **solver_kwargs)

    return -problem.value, problem.solver_stats


def qubo_sdp_via_maxcut_01(Q, solver=cvx.MOSEK, verbose=False, **solver_kwargs):
    """
    Compute sdp relaxation of QUBO via MAX-CUT.
    solvers = [cvx.SCS, cvx.CVXOPT, cvx.MOSEK]
    MOSEK licence: https://www.mosek.com/products/version-11/
    """
    MC = qubo2maxcut(Q) # Convert to a MAX-CUT 

    return maxcut_sdp(MC, solver=solver, verbose=verbose, **solver_kwargs)


def qbo_sdp_01(Q, linear_constraint=None, quadratic_constraints=None, 
               solver=cvx.MOSEK, verbose=False, **solver_kwargs):
    """Compute SDP relaxation of QBO problem x^TQx, x in {0, 1}^n

    Args:
        Q (square numpy 2d array): A matrix Q.
        linear_constraint (Tuple(A, b, operator), optional): A linear constraint. Defaults to None.
        quadratic_constraints (Iterable(Tuple(Qi, ri)), optional): A collection of quadratic constraints. Defaults to None.
        solver (cvx.solver, optional): A solver for the SDP relaxation. Defaults to cvx.MOSEK.
        verbose (bool, optional): Verbosity of cvx.solver. Defaults to False.

    Raises:
        TBD
    Returns:
        TBD
    """
    # Allowed (in)equality operators in constraints
    allowed_operators = {'<=': operator.le, '>=': operator.ge, '==': operator.eq}

    def check_operator(operator_as_string):
        """Check operator.

        Args:
            operator_as_string (string): A string representation of (in)equality operator - python syntax.

        Raises:
            ValueError: Raises ValueError of not allowed operator is used.

        Returns:
            TBD
        """
        try:
            return allowed_operators[operator_as_string]
        except KeyError:
            raise ValueError(f'Operator {operator_as_string} is not allowed. Please use one of {list(allowed_operators.keys())}')
        except Exception as ex:
            raise(ex)
        
    
    def check_linear_constraint(linear_constraint, n):
        """TBD

        Args:
            linear_constraint (_type_): _description_
            n (_type_): _description_

        Raises:
            ValueError: _description_
            ex: _description_
            ValueError: _description_
            ex: _description_
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        A, b, operator_as_string = linear_constraint
        _operator = check_operator(operator_as_string)
        try:
            m_A, n_A = A.shape
        except ValueError:
            raise ValueError(f'A has to be a numpy 2D array of a shape (m, {n}). Got: {A.shape}')
        except Exception as ex:
            raise ex()
            
        try:
            m_b, = b.shape
        except ValueError:
            raise ValueError(f'b has to be a numpy 1D array. Got: {b.shape}')
        except Exception as ex:
            raise ex()

        if not (n == n_A or m_A == m_b):
            raise ValueError(f'Incompatible dimensions. n={n}, A=({m_A}, {n_A}), b=({m_b})')

        return A, b, _operator

    def check_quadratic_constraints(quadratic_constrains, n):
        """TBD

        Args:
            quadratic_constrains (_type_): _description_
            n (_type_): _description_

        Raises:
            ValueError: _description_
            ex: _description_
            ValueError: _description_

        Yields:
            _type_: _description_
        """
        for i, (Qi, ri, operator_as_string) in enumerate(quadratic_constraints):
            _operator = check_operator(operator_as_string)
            try:
                m_Qi, n_Qi = Qi.shape
            except ValueError:
                raise ValueError(f'Q{i} has to be a square numpy 2D array of a shape ({n}, {n}). Got: {Qi.shape}')
            except Exception as ex:
                raise ex()   

            if not (m_Qi == n_Qi and n == m_Qi):
                raise ValueError(f'Q{i} has to be a square numpy 2D array of a shape ({n}, {n}). Got: {Qi.shape}')

            yield Qi, ri, _operator
            
    
    start_time = time.time()

    n, _ = Q.shape
    
    Y = cvx.Variable((n+1, n+1), symmetric=True)
    X = Y[1:, 1:]
    x = cvx.diag(X)
    constraints = [
        Y >> 0, 
        Y[0, 0]== 1, 
        Y[0, 1:] == x, 
        x >= 0, x <= 1
    ]

    objective = cvx.Minimize(cvx.trace(Q @ X))

    if linear_constraint is not None:
        A, b, _operator = check_linear_constraint(linear_constraint, n)
        A_m, _ = A.shape
        
        constraints += [
            _operator(A[i, :] @ x, b[i]) for i in range(A_m)
        ]


    if quadratic_constraints is not None:
        constraints += [
            _operator(cvx.trace(Qi @ X), ri) for Qi, ri, _operator in check_quadratic_constraints(quadratic_constraints, n)
        ]
        
    problem = cvx.Problem(objective, constraints)
    problem.solve(solver=solver, verbose=verbose, **solver_kwargs)
    runtime = time.time() - start_time

    
    return (problem.value,
            (x.value, X.value, Y.value, problem.status, runtime))


def transform_from_01(Y):
    """
    Transforms the matrix Y, which is in 0/1 encoding, to the matrix Z, which is in -1/1 encoding.
    This can then be used for the Goemans-Williamnson rounding routine. 
    n is the problem size
    """
    N, _ = Y.shape # is of size n + 1
    n = N - 1
    
    # Step 1: Construct vector e of size n = N - 1
    e = np.ones(n)
    
    # Step 2: Construct the matrix V (size (n+1) x (n+1))
    I_n = np.eye(n)  # Identity matrix of size n x n
    V = np.zeros((N, N))
    
    # Fill in the parts of the matrix V
    V[0, 0] = 1           # Top-left block is 1
    V[1:, 0] = -e         # First column (without Y[0,0]) is -e
    V[1:, 1:] = 2*I_n     # Top-right block is zero
    
    # Step 3: Perform matrix multiplication to get Z
    Z = V @ Y @ V.T
    
    return Z


def transform_to_01(Z):
    """
    Transforms the matrix Z, which is in -1/1 encoding, to the matrix Y, which is in 0/1 encoding.
    This can then be then used to extract the results from the Goemans-Williamnson rounding routine. 
    n is the problem size
    """
    N, _ = Z.shape
    n = N - 1

    e = np.ones(n)
    I_n = np.eye(n)
    U = np.zeros((N, N))

    U[0, 0] = 1
    U[1:, 0] = 0.5*e
    U[1:, 1:] = 0.5*I_n

    Y = U @ Z @ U.T

    return Y


def gw_unit_vectors_from_Z(Z):
    """
    Factor a PSD matrix Z (diag = e) into rows of unit vectors
    Do an eigenvalue decomposition, rebuild V = U sqrt(Λ) and row-normalize for
    numerical robustness
    """
    #Z = 0.5 * (Z + Z.T)               # symmetrize - not needed, matrix is symmetric due to definition of the problem
    w, U = np.linalg.eigh(Z)          # Z = U Diag(w) U^T = U Diag(sqrt(w)) Diag(sqrt(w)) V^T
    w = np.clip(w, 0.0, None)         # Clip tiny negative eigenvalues and map them to zero
    V = U @ np.diag(np.sqrt(w))       # V V^T ~ Z

    # Row-normalize so rows are ~ unit length
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return V / norms


def gw_round_once(V, rng=None):
    """
    One GW hyperplane rounding from Z (+-1 correlation SDP in lifted space)
    Returns a +-1 vector s of length N (index 0 is the anchor)
    V is the decomposition of Z, we do that only once
    """
    if rng is None:
        rng = np.random.default_rng()
    n = V.shape[1]
    r = rng.normal(size=n) # random hyperplane
    #r = np.random.randn(n)
    r = r/np.linalg.norm(r) # normalized, from the unit sphere
    #print("Random hyperplane: ", r)
    s = np.sign(V @ r)
    #print("Signs: ", s)
    s[s == 0] = 1                     
    return s


def max_viol_linear(x, linear_constraint):
    """
    Returns the maximum violation of linear constraints
    """
    if linear_constraint is None:
        return 0.0
    else:
        A, b, _op = linear_constraint
        viol = np.linalg.norm(A @ x - b, np.inf)
        return viol


def max_viol_quadratic(x, constraint_list=None):
    """
    constraints_list is an iterable of (Qi, ri, op_string). We report
    max absolute violation x^T Qi x - ri (operator-agnostic).
    """
    if constraint_list is None: 
        return 0.0 # No violation if no constraints exist
    else:
        v = 0.0
        for (Qi, ri, _op) in constraint_list:
            #print(f"Q_i, r_i, _op: {Qi, ri, _op}")
            Qi = np.asarray(Qi)
            # Ensure x is evaluated numerically (after solving the optimization)
            v = max(v, abs(float(np.vdot(x.T, (Qi @ x)) - ri)))
        return v


def goemans_williamson_loop(
    Y, Q, linear_constraint = None, quadratic_constraints = None, tol_lin = 1e-6, tol_quad = 1e-6,
    max_trials = 100, seed = None,
    local_improve = None,
    stop_at_first_feasible = False,
):
    """
    GW rounding loop in the BiqCrunch spirit (minimization setting).

    Inputs
    ------
    Y : (n+1)x(n+1) PSD 'lifted' matrix in {0,1}-encoding 
    Q : (n)x(n) symmetric objective matrix; objective is x^T Q x (to be MINIMIZED)
    linear_constraint : optional triple (A, b, meta) representing A x <= b (or similar)
    quadratic_constraints : optional structure consumed by max_viol_quadratic(x, ...)
    tol_lin : feasibility tolerance for linear constraints
    tol_quad : feasibility tolerance for quadratic constraints
    max_trials : number of random hyperplanes to try
    seed : RNG seed
    local_improve : optional callable that maps a feasible x to (x_loc, val_loc) with
                    val_loc <= x^T Q x and x_loc feasible; run until locally optimal at the end of the GW computation
    stop_at_first_feasible : if True, return the first feasible x (after local improve if provided)

    Returns
    -------
    (x_star, info)
        x_star : best feasible {0,1}^n found (None if none found).
        info   : diagnostics dict with keys:
                 {
                   'accepted': bool,
                   'trials_used': int,
                   'ub_value': float or None,
                   'best_infeasible': {'value','x','vlin','vquad','trial'} or None
                 }
    """
    rng = np.random.default_rng(seed)

    # 0/1-lifted Y -> {-1,1}-lifted Z; decompose once to get unit vectors for GW
    Z = transform_from_01(np.asarray(Y))                  # shape (n+1, n+1)
    V = gw_unit_vectors_from_Z(Z)                         # rows are unit vectors; shape (n+1, d)

    # Unpack linear constraint, if any
    if linear_constraint is not None:
        A, b, _op = linear_constraint
        A = np.asarray(A)
        b = np.asarray(b)

    n = Q.shape[0]

    x_star = None            # incumbent feasible
    beta = np.inf            # incumbent UB value (since we're minimizing)
    best_infeasible = None   # diagnostics only

    for t in range(1, max_trials + 1):
        # One GW round: produce s in {-1,1}^{n+1}, then map to x in {0,1}^n
        s = gw_round_once(V, rng)  # length n+1, s[0] is the anchor
        x = (s[1:] + 1.0) * 0.5     # length n, in {0,1}

        # Violations
        vlin = 0.0 if linear_constraint is None else max_viol_linear(x, linear_constraint)
        vquad = max_viol_quadratic(x, quadratic_constraints)

        feasible = (vlin <= tol_lin) and (vquad <= tol_quad)

        # Objective (upper bound)
        val = float(np.vdot(x, Q @ x))  # x^T Q x

        if feasible:
            # Optional local improvement in the feasible region, BiqCrunch-style
            # TODO: move this to the end of the loop, this is only done once after obtaining a feasible solution for the original problem
            if local_improve is not None:
                x_loc, val_loc = local_improve(x)
                # Trust but verify feasibility after local_improve
                vlin_loc = 0.0 if linear_constraint is None else max_viol_linear(x_loc, linear_constraint)
                vquad_loc = max_viol_quadratic(x_loc, quadratic_constraints)
                if (vlin_loc <= tol_lin) and (vquad_loc <= tol_quad) and (val_loc <= val):
                    x, val, vlin, vquad = x_loc, val_loc, vlin_loc, vquad_loc

            # BiqCrunch "test of improvement": accept only if it improves UB
            if val < beta:  # strict improvement
                x_star, beta = x, val

            if stop_at_first_feasible:
                return x_star, {
                    'accepted': True,
                    'trials_used': t,
                    'ub_value': beta,
                    'best_infeasible': best_infeasible,
                    'vlin': vlin,
                    'vquad': vquad
                }
        else:
            # Track best infeasible for reporting
            if best_infeasible is None or val < best_infeasible['value']:
                best_infeasible = {
                    'value': val,
                    'x': x.copy(),
                    'vlin': vlin,
                    'vquad': vquad,
                    'trial': t
                }

    # End of trials: Always return the best found solution (feasible or infeasible)
    if x_star is not None:
        return x_star, {
            'accepted': True,
            'trials_used': max_trials,
            'ub_value': beta,
            'best_infeasible': best_infeasible,
            'vlin': 0.0,  # no violations if feasible
            'vquad': 0.0
        }
    else:
        # No feasible solution found: return best infeasible
        return None, {
            'accepted': False,
            'trials_used': max_trials,
            'ub_value': beta,
            'best_infeasible': best_infeasible,
            'vlin': best_infeasible['vlin'] if best_infeasible else np.nan,
            'vquad': best_infeasible['vquad'] if best_infeasible else np.nan
        }


def load_data(filename):
    with open(filename) as f:
        data = json.load(f)


    Q = np.asarray(data['QBO']['Q'])
    linear_constraint = data['QBO']['constraints']['linear']
    quadratic_constraints = data['QBO']['constraints']['quadratic']
    if linear_constraint is not None:
        A, b, _operator = linear_constraint
        linear_constraint = np.asarray(A), np.asarray(b), _operator

    if quadratic_constraints:
        quadratic_constraints = ((np.asarray(Qi), ri, _operator) for Qi, ri, _operator in quadratic_constraints)

    return Q, linear_constraint, quadratic_constraints, data

if __name__ == "__main__":

    filename = sys.argv[1]
    print('Computing problem instance', filename)

    Q, linear_constraint, quadratic_constraints, data = load_data(filename)

    (value, (x, X, Y, status, runtime)) = qbo_sdp_01(Q, 
                                                     solver=cvx.CVXOPT, 
                                                     linear_constraint=linear_constraint, 
                                                     quadratic_constraints=quadratic_constraints, 
                                                     verbose=True)
    
    gw_x, gw_info = goemans_williamson_loop(Y, Q,
                                            linear_constraint=linear_constraint,
                                            quadratic_constraints=quadratic_constraints,
                                            tol_lin=1e-6, tol_quad=1e-6,
                                            max_trials=10, seed=42,
                                            local_improve=None,               # or "hill-climber", function for local improvement, TBD
                                            stop_at_first_feasible=False
)

    print('Computed lower bound: ', value)
    print('Optimal value:', data['optimum'])

    print("\n[GW] trials used:", gw_info["trials_used"])

    if gw_info["accepted"] and gw_x is not None:
        # Recompute violations for the returned incumbent (safe + explicit)
        vlin = 0.0 if linear_constraint is None else max_viol_linear(gw_x, linear_constraint)
        vquad = max_viol_quadratic(gw_x, quadratic_constraints)

        print("[GW] feasible solution found")
        print("[GW] upper bound (x^T Q x):", gw_info["ub_value"])
        print("[GW] max linear violation:", vlin)
        print("[GW] max quadratic violation:", vquad)
        print("[GW] x (0/1):", gw_x)
    else:
        print("[GW] no feasible solution found")
        bi = gw_info.get("best_infeasible", None)
        if bi is not None:
            print("[GW] best infeasible (diagnostic) at trial:", bi["trial"])
            print("[GW] objective (x^T Q x):", bi["value"])
            print("[GW] max linear violation:", bi["vlin"])
            print("[GW] max quadratic violation:", bi["vquad"])
            print("[GW] x (0/1):", bi["x"])


