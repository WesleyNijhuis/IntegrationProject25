function u0 = Sget_mu(x,reference,t)

    persistent fast_opt

    if t==0
        dt=0.01;
        mpc_A = [0.996,0.005906,0.00976,-0.00005518;0.000003775,1.012,-0.0001249,0.009988;-0.08811,1.168,0.9523,-0.008953;0.001128,2.431,-0.02479,1.002];
        mpc_B = [0.1022;0.0961;20.25;19.14];
        mpc_C = [1,0,0,0;0,1,0,0];
        mpc_D = [0;0];
        
        Q_mpc = diag([1e1, 1e-1, 1e-3,1e-3]);
        R_mpc = 1e-1;

        yalmip('clear')
        
        [P,K,~] = idare(mpc_A,mpc_B,Q_mpc,R_mpc);
        
        % mptopt('qpsolver', 'quadprog');
        % mptopt('lpsolver','LCP')
        
        mpt3_model = ss(mpc_A,mpc_B,mpc_C,mpc_D,dt);
        mpc_mpt3 = LTISystem(mpt3_model);
        mpc_mpt3.x.min = [-pi/2; -pi/4; -inf; -inf];
        mpc_mpt3.x.max = [pi/2; pi/4; inf; inf];
        mpc_mpt3.u.min = -1;
        mpc_mpt3.u.max = 1;
        mpc_mpt3.x.with('reference');
        mpc_mpt3.x.reference = 'free';
        
        mpc_mpt3.x.penalty = QuadFunction(Q_mpc);
        mpc_mpt3.u.penalty = QuadFunction(R_mpc);
        Tset = mpc_mpt3.LQRSet();
        PN = mpc_mpt3.LQRPenalty;
        mpc_mpt3.x.with('terminalSet');
        mpc_mpt3.x.terminalSet = Tset;
        mpc_mpt3.x.with('terminalPenalty');
        mpc_mpt3.x.terminalPenalty = QuadFunction(P);
        
        mpt_horizon = 20;
        ctrl = MPCController(mpc_mpt3,mpt_horizon);%.toExplicit();    
    
        
        % Export to a (faster) MPC controller
        %ctrl.optimizer.toMatlab('mycontroller.m', 'primal', 'obj');
        
        % precompile optimizer and using YALMIP directly for speed
        optim = ctrl.toYALMIP();
        constr = optim.constraints;
        obj    = optim.objective;
        vars   = optim.variables;
        params = [vars.x(:,1); vars.filters.x.reference];
        decision = vars.u;
        fast_opt = optimizer(constr, obj, sdpsettings('solver','quadprog'), params, decision);
    
        mu = fast_opt{[x; reference]};
        u0 = mu(1)
        
        if u0 == NaN
            u0 =0;
        end
    else
        fast_opt = optimizer(constr, obj, sdpsettings('solver','quadprog'), params, decision);
    
        mu = fast_opt{[x; reference]};
        u0 = mu(1);
        
        if u0 == NaN
            u0 =0;
        end
    end


end