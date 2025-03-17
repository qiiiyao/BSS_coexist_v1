def compute_scheme_for_LV_focal(A, focal_sp, mu = None):
    """
    Computes the invasion scheme for a LV community model
    
    Computes all feasible subcommunities and the corresponding invasion growth
    rates for a community model of the form
        1/N_i dN_i/dt = \mu_i - \sum_j a_ij N_j
    

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates

    Returns
    -------
    inv_scheme : numpy array, shape (n_equi, n_spec)
        Invasion scheme containing the invasion growth rates of the ``n_spec``
        species in the ``n_equi`` different equilibrium densities
        Non-stable equilibria may be included with np.nan rows
    """
    A = pd.read_csv("results/fit_results/BSS_exculde_trees_raw/parameters/a_10_1_c.csv", index_col = 0).values
    n = len(A)
    focal_sp = 3
    if mu is None:
      mu = np.ones(n)
    counter = 0
    equi_all = np.zeros(( 2**n, n))
    r_i = np.full(equi_all.shape, np.nan)
    # compute all possible equilibria
    for n_spec in range(1, n+1):
      combs = np.array(list(combinations(range(n), n_spec)))
      combs = [arr for arr in combs if focal_sp in arr]

      A_temp = np.empty((len(combs), n_spec, n_spec))
      mu_temp = np.empty((len(combs), n_spec))
      for i, comb in enumerate(combs):
        A_temp[i] = A[comb[:,np.newaxis], comb]
        mu_temp[i] = mu[comb]
        
        # solve equilibrium vectorized for speed
      equi_temp = np.linalg.solve(A_temp, mu_temp)
        
      for i, comb in enumerate(combs):
        counter += 1
        equi_all[counter, comb] = equi_temp[i]
  # community when all species are absent
    equi_all[0] = 0
     
  # delete all zero row
    zero_rows = np.all(equi_all == 0, axis=1)
    equi_all = equi_all[~zero_rows]
  # compute invasion growth rates
    r_i = mu - np.sum(A*equi_all[:,np.newaxis], axis = -1)
    r_i = r_i[~np.any(equi_all<0, axis = 1)]
    
  # remove rounding errors
    r_i_rounded = np.round(r_i, 10)
    return r_i_rounded, equi_all

def equi_for_invasion_focal(A, focal_sp, mu = None, remove_non_feasible = True,
                      ret_groups = False):
    """For each species computes the resident community at the end state
    
    Given the LV-community  given by
        1/N_i dN_i/dt = \mu_i - \sum_j a_ij N_j
    it computes the invasion graph and a possible end state of assembly.
    Given this end state of assembly it computes the resident community for all
    species at invasion.
    For the species not present at the end-state, the resident community is
    simply the end state of assembly.
    For the species present the resident community is the n-1 community of
    the end state of assembly. If this n-1 community does not exist then the
    resident community is np.nan.

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates. Default is mu = 1
    remove_non_feasible : Bool, optional
        Whether non-feasible resident communities are replaced with np.nan.
        Recommended unless in debugging. The default is True.
    ret_groups : bool, optional
        Whether to return the different species groups along with the resident
        communities. If ''True'', then excl, w_n1 and wo_n1 are returned.
        The default is False.

    Raises
    ------
    RuntimeError
        If there are more than one end state to assembly

    Returns
    -------
    equi_invasion: numpy array, shape (n_spec, n_spec)
        equi_invaison[i] is the density of the resident community when species
        i invades
    excl: numpy array
        indices of species which are excluded from the end state of assembly.
        Only returned if ret_groups is True
    wo_n1: numpy array:
        Indices of species which are present at the end state of assembly, but
        do not have an n-1 community. For these species equi_invasion[i] is np.nan
        Only returned if ret_groups is True
    wo_n1: numpy array:
        Indices of species which are present at the end state of assembly, and
        do have an n-1 community. For these species equi_invasion[i] is np.nan
        Only returned if ret_groups is True

    """
    if mu is None:
        mu = np.ones(len(A))
    #focal_sp = 3
    # compute equilibrium and invasion growth rates
    A = pd.read_csv("results/fit_results/BSS_exculde_trees_raw/parameters/a_10_1_c.csv", index_col = 0).values
    r_i, equi_all = compute_scheme_for_LV(A)
    adj_real, cyclic, permanent = compute_graph(r_i)
    if len(permanent) != 1:
        raise RuntimeError("Feature not implmented yet")
    equi_all_c = [arr for arr in equi_all if all(item >= 0 for item in arr)]
    specific_position = 3  # ????????????
    min_zero_count = np.inf  # ??????????????????0??????????????????0?????????
    selected_array = None  # ???????????????????????????

    # ????????????????????????
    for arr in equi_all_c:
    # ??????0?????????
        zero_count = np.count_nonzero(arr == 0)
    
    # ???????????????????????????????????????0
        if zero_count < min_zero_count and arr[specific_position-1] == 0:
            min_zero_count = zero_count
            selected_array = arr
    non_zero_indices = np.nonzero(selected_array) 

    # species present
    pres = np.where(r_i[permanent[0]] == 0)[0]
    # species which are excluded
    excl = np.where(r_i[permanent[0]] != 0)[0]
    
    # compute the resident community for all invasions
    equi_invasion = np.zeros((len(mu), len(mu)))
    equi_invasion[excl[:,np.newaxis], pres] = np.linalg.solve(A[pres[:,np.newaxis], pres], mu[pres])
    
    # compute invasion densities for the other species
    for i in pres:
        i = 3
        sub_p = np.delete(non_zero_indices, i-1)
        sub_p = np.concatenate(non_zero_indices, axis=0)
        equi_invasion[i, sub_p] = np.linalg.solve(A[sub_p[:,np.newaxis], sub_p],
                                             mu[sub_p])
    
    if remove_non_feasible:
        # some equilibria might not be feasible
        equi_invasion[np.any(equi_invasion<0, axis = 1)] = np.nan
    
    if ret_groups:
        return equi_invasion, excl, pres[np.isnan(equi_invasion[pres,0])], pres[np.isfinite(equi_invasion[pres,0])]
    return equi_invasion
