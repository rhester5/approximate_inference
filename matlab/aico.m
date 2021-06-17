% clear all
close all

pendulum = true;
cartpole = false;
drones = false;

if pendulum
    xn = 2;
    un = 1;
elseif cartpole
    xn = 4;
    un = 1;
elseif drones
    xn = 12;
    un = 6;
end

b = 0.3;
k = 1;
l = 1;
num_iter = 10;

Q = eye(xn); %*100; % kind of works with * 1
R = eye(un);
Qf = eye(xn); %*100; % kind of works * 1

q = eye(xn);
r = eye(un);

num_steps = 100;
dt = 0.1;

% compute desired trajectory using direct collocation

if pendulum
    x0 = [0, 0];
    u0 = [0];
    xd = [pi, 0];
    xd2 = -xd;
elseif cartpole
    x0 = [0, 0, 0, 0];
    u0 = [0];
    xd = [0, pi, 0, 0];
    xd2 = -xd;
elseif drones
    % x0 = [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
    x0 = [-2, -2, 1, 1, 4, 4, 0, 0, 0, 0, 0, 0];
    u0 = zeros(1, un);
    xd = [0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
end

% [x_d, u_d] = direct_collocation(x0, xd, u0, q, r, num_steps, dt, b, k, l, pendulum, cartpole, drones);
% x_d = get_xd(x0, xd, dt, num_steps);
u_d = zeros(un, 100);

% TODO why aren't I using x_d anywhere
% TODO because you need to add the cost function shit

disp('Computed desired trajectory');

% compute nominal trajectory using direct collocation

% [x0, u0] = direct_collocation(x0, xd2, u0, q, r, num_steps, dt, b, k, l, pendulum, cartpole, drones);
x0 = zeros(xn, 100);
u0 = zeros(un, 100);

disp('Computed nominal trajectory');

x = cell([1, num_iter+1]);
u = cell([1, num_iter+1]);

for it=1:num_iter+1
   x{it} = x0;
   u{it} = u0; 
end

q = diag(Q);
r = diag(R);
tol = 1e-3;
thresh = 1e-3;
max_step_size = 0.1;
mu = zeros(xn, 1);
Cov = eye(xn) * 1e0;
advance_b_beyond_xhat = true;
damping = b
Damping = 1;
max_reloc_fwd = 0;
break_when = true;
Damping_ref = cell([1, num_steps+1]); % TODO gonna roll with Damping_ref being num_steps not num_steps+1 (but for now changing it in update belief)
for i=1:num_steps+1
    Damping_ref{i} = 1;
end

sweep = 0;
scale = 0;
c = -1;
useBwdMsg = false;
fixFinalState = false;

% init messages
% messages
  % s.resize(T+1, n);  Sinv.resize(T+1, n, n);
  % v.resize(T+1, n);  Vinv.resize(T+1, n, n);
  % r.resize(T+1, n);  R.resize(T+1, n, n);     r.setZero();  R   .setZero();
  % for(uint t=0;t<=T;t++){
  %   s[t]=x0;  Sinv[t].setDiag(1e-10);
  %   v[t]=x0;  Vinv[t].setDiag(1e-10);
  % }
  % s[0]=x0;  Sinv[0].setDiag(1e10);
  % if(useBwdMsg){ v[T] = bwdMsg_v;  Vinv[T] = bwdMsg_Vinv; }

s = cell([1, num_steps+1]);
Sinv = cell([1, num_steps+1]);
v = cell([1, num_steps+1]);
Vinv = cell([1, num_steps+1]);
% r = cell([1, num_steps+1]);
% R = cell([1, num_steps+1]);
for i=1:num_steps+1
    s{i} = x0(:, 1);
    Sinv{i} = eye(xn) * 1e-10;
    v{i} = x0(:, 1);
    Vinv{i} = eye(xn) * 1e-10;
    % r{i} = zeros(xn, 1);
    % R{i} = zeros(xn, xn);
end
Sinv{1} = eye(xn) * 1e10;
if useBwdMsg
    v{num_steps+1} = bwdMsg_v;
    Vinv{num_steps+1} = bwdMsg_Vinv;
    % TODO should this be num_steps or num_steps+1
end

% beliefs
  % b.resize(T+1, n);  Binv.resize(T+1, n, n);  b.setZero();  Binv.setZero();  b[0]=x0;  Binv[0].setDiag(1e10);
  % rhat.resize(T+1);     rhat.setZero();
  % xhat.resize(T+1, n);  xhat.setZero();  xhat[0]=x0;
  % //if(!sys->isKinematic()) soc::getPositionTrajectory(q, b); else q=b;
  % //dampingReference = qhat;
  % dampingReference.clear();

b = cell([1, num_steps+1]);
Binv = cell([1, num_steps+1]);
% rhat = cell([1, num_steps+1]);
xhat = cell([1, num_steps+1]);
for i=1:num_steps % TODO should this be num_steps+1 or should the sizes be num_steps
    b{i} = zeros(xn, 1);
    % b{i} = x_d(:, i);
    Binv{i} = zeros(xn, xn);
    % rhat{i} = zeros(1, 1);
    % xhat{i} = zeros(xn, 1);
    % xhat{i} = x0(:, i);
    xhat{i} = x_d(:, i);
end
xhat{1} = x0(:, 1);

% dynamics aren't affine so no need for state dynamics vector
a = zeros(xn, 1);

tStart = cputime;
tIter = tStart;

for it = 1:num_iter

msg = ['AICO iteration #', num2str(it)];
disp(msg);

% rememberOldState
b_old = b;
xhat_old = xhat;
s_old = s;
Sinv_old = Sinv;
v_old = v;
Vinv_old = Vinv;

for t = 1:num_steps-1
    % update Time Step
    if pendulum
        [A, B] = linear_pendulum_dynamics(x{it}(:, t), damping);
    elseif cartpole
        A = A_cartpole(x{it}(:, t), u{it}(:, t));
        B = B_cartpole(x{it}(:, t), u{it}(:, t));
    end

    if t > 1
        % update Fwd Message
        % [A, B] = linear_pendulum_dynamics(x{it}(:, t-1), damping);
        % arr barS, St;
        % inverse_SymPosDef(barS, Sinv[t-1] + R[t-1]);
        barS = inv(Sinv{t-1} + Q);
        % St = Q[t-1];
        St = Cov;
        % St += B[t-1]*Hinv[t-1]*tB[t-1];
        St = St + B * inv(R) * B';
        % St += A[t-1]*barS*tA[t-1];
        St = St + A * barS * A';
        % s[t] = a[t-1] + A[t-1]*(barS*(Sinv[t-1]*s[t-1] + r[t-1]));
        % a = mvnrnd(mu, Cov);
        % a = reshape(a, xn, 1);
        s{t} = a + A * (barS * (Sinv{t-1} * s{t-1} + q));
        % inverse_SymPosDef(Sinv[t](), St);
        Sinv{t} = St;
    end

    % update Belief
    % if(damping && dampingReference.N){
    % TODO is this what dampingReference.N means
    if Damping && norm(Damping_ref{num_steps})
    %     Binv[t] = Sinv[t] + Vinv[t] + R[t] + damping*eye(R.d1);
        Binv{t} = Sinv{t} + Vinv{t} + Q + Damping * eye(xn);
    %     lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t] + damping*dampingReference[t]);
        b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q + Damping * Damping_ref{t});
        % TODO is this actually what Ainv b is implying
    %   }
    else
        % Binv[t] = Sinv[t] + Vinv[t] + R[t];
        Binv{t} = Sinv{t} + Vinv{t} + Q;
        % lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t]);
        b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q);
    end

    if it > 1
        max_reloc = max_reloc_fwd;
    else
        max_reloc = 0;
    end
    for k = 1:max_reloc
        % if(k || !forceRelocation) where forceRelocation is true if cost < 0
        if k > 0 || c > 0
        %   if(maxDiff(b[t], xhat[t])<tolerance) break;
            if abs(b{t} - xhat{t}) < tol
                break
            end
        end

        % update Task Message
        % xhat_t = b[t]
        xhat_t = b{t};
        % if(maxStepSize>0. && norm(xhat_t-xhat[t])>maxStepSize){
        if norm(xhat_t - xhat{t}) > max_step_size
        %     arr Delta = xhat_t-xhat[t];
            delta = xhat_t - xhat{t};
        %     Delta *= maxStepSize/norm(Delta);
            delta = delta * max_step_size / norm(delta);
        %     xhat_t = xhat[t] + Delta;  //really change the given xhat_t (often the belief!!)
            xhat_t = xhat{t} + delta;
        % }
        end
        % xhat[t]() = xhat_t;
        xhat{t} = xhat_t;
        % countSetq++;
        % sys->setx(xhat[t]);
        x{it+1}(:, t) = xhat{t};

        if advance_b_beyond_xhat
            % update Belief
            % if(damping && dampingReference.N){
            % TODO is this what dampingReference.N means
            if Damping && norm(Damping_ref{num_steps})
            %     Binv[t] = Sinv[t] + Vinv[t] + R[t] + damping*eye(R.d1);
                Binv{t} = Sinv{t} + Vinv{t} + Q + Damping * eye(xn);
            %     lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t] + damping*dampingReference[t]);
                b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q + Damping * Damping_ref{t});
                % TODO is this actually what Ainv b is implying
            %   }
            else
                % Binv[t] = Sinv[t] + Vinv[t] + R[t];
                Binv{t} = Sinv{t} + Vinv{t} + Q;
                % lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t]);
                b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q);
            end
        end
    end
end

disp('Completed forward pass');

for i = 1:num_steps
    t = num_steps-i+1;
    % update Time Step

    if pendulum
        [A, B] = linear_pendulum_dynamics(x{it}(:, t), damping);
    elseif cartpole
        A = A_cartpole(x{it}(:, t), u{it}(:, t));
        B = B_cartpole(x{it}(:, t), u{it}(:, t));
    end

    if ~(fixFinalState && t == num_steps)
        if t < num_steps
            % update Bwd Message
            % arr barV, Vt;
            % inverse_SymPosDef(barV, Vinv[t+1] + R[t+1]);
            barV = inv(Vinv{t+1} + Q);
            % Vt = Q[t];
            Vt = Cov;
            % Vt += B[t]*Hinv[t]*tB[t];
            Vt = Vt + B * inv(R) * B';
            % Vt += barV;
            Vt = Vt + barV;
            % Vt = Ainv[t]*Vt*invtA[t];
            Vt = inv(A) * Vt * inv(A');
            % v[t] = Ainv[t]*(-a[t] + barV*(Vinv[t+1]*v[t+1] + r[t+1]));
            % a = mvnrnd(mu, Cov);
            % a = reshape(a, xn, 1);
            v{t} = inv(A) * (-a + barV * (Vinv{t+1} * v{t+1} + q));
            % inverse_SymPosDef(Vinv[t](), Vt);
            Vinv{t} = inv(Vt);
        else
            if ~useBwdMsg
                % v[t] = b[t];
                v{t} = b{t};
                % Vinv[t].setDiag(1e-4); 
                Vinv{t} = Vinv{t} - diag(diag(Vinv{t})) + eye(xn) * 1e-4;
            else
                % v[T] = bwdMsg_v;
                v{num_steps+1} = bwdMsg_v;
                % TODO num_steps or num_steps+1
                % Vinv[T] = bwdMsg_Vinv;
                Vinv{num_steps+1} = bwdMsg_Vinv;
            end
        end
    end

    % update Belief
    % if(damping && dampingReference.N){
    % TODO is this what dampingReference.N means
    if Damping && norm(Damping_ref{num_steps})
    %     Binv[t] = Sinv[t] + Vinv[t] + R[t] + damping*eye(R.d1);
        Binv{t} = Sinv{t} + Vinv{t} + Q + Damping * eye(xn);
    %     lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t] + damping*dampingReference[t]);
        b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q + Damping * Damping_ref{t});
        % TODO is this actually what Ainv b is implying
    %   }
    else
        % Binv[t] = Sinv[t] + Vinv[t] + R[t];
        Binv{t} = Sinv{t} + Vinv{t} + Q;
        % lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t]);
        b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q);
    end

    max_reloc = 1;
    for k = 1:max_reloc
        % if(k || !forceRelocation) where forceRelocation is true if cost < 0
        if k > 0 || c > 0 % TODO do I need to use a different variable name for cost
        %   if(maxDiff(b[t], xhat[t])<tolerance) break;
            if abs(b{t} - xhat{t}) < tol
                break
            end
        end

        % update Task Message
        % xhat_t = b[t]
        xhat_t = b{t};
        % if(maxStepSize>0. && norm(xhat_t-xhat[t])>maxStepSize){
        if norm(xhat_t - xhat{t}) > max_step_size
        %     arr Delta = xhat_t-xhat[t];
            delta = xhat_t - xhat{t};
        %     Delta *= maxStepSize/norm(Delta);
            delta = delta * max_step_size / norm(delta);
        %     xhat_t = xhat[t] + Delta;  //really change the given xhat_t (often the belief!!)
            xhat_t = xhat{t} + delta;
        % }
        end
        % xhat[t]() = xhat_t;
        xhat{t} = xhat_t;
        % countSetq++;
        % sys->setx(xhat[t]);
        x{it+1}(:, t) = xhat{t};

        if advance_b_beyond_xhat
            % update Belief
            % if(damping && dampingReference.N){
            % TODO is this what dampingReference.N means
            if Damping && norm(Damping_ref{num_steps})
            %     Binv[t] = Sinv[t] + Vinv[t] + R[t] + damping*eye(R.d1);
                Binv{t} = Sinv{t} + Vinv{t} + Q + Damping * eye(xn);
            %     lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t] + damping*dampingReference[t]);
                b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q + Damping * Damping_ref{t});
                % TODO is this actually what Ainv b is implying
            %   }
            else
                % Binv[t] = Sinv[t] + Vinv[t] + R[t];
                Binv{t} = Sinv{t} + Vinv{t} + Q;
                % lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t]);
                b{t} = inv(Binv{t}) * (Sinv{t} * s{t} + Vinv{t} * v{t} + q);
            end
        end

    end
end

disp('Completed backward pass');

% b_step=maxDiff(b_old, b);
% dampingReference=b;
Damping_ref = b;

% if(!advanceBeliefBeyondXhat || sweepMode==smILQG){
% cost = evaluateTrajectory(b, display>0); //this routine takes the current R, r matrices to compute costs
% % double exact_cost = sys->analyzeTrajectory(b, display>0); //this routine calles the simulator again for each time step
% % sys->costChecks(b);
% % cout <<"DIFF=" <<fabs(cost-exact_cost) <<endl;
% }else{
% cost = analyzeTrajectory(*sys, b, display>0, &cout); //this routine calles the simulator again for each time step
% }
c = cost(x{it+1}, u{it+1}, x_d, u_d, Q, R, Qf, num_steps);

% if(damping) perhapsUndoStep();
if Damping
% if(cost_old>0 && cost>cost_old){
    if c > cost(x{it}, u{it}, x_d, u_d, Q, R, Qf, num_steps)
%     //cout <<"\b AICO REJECT: cost=" <<cost <<" cost_old=" <<cost_old <<endl;
%     damping *= 10.;
    Damping = Damping * 10;
%     dampingReference = b_old;
    Damping_ref = b_old;
%     cost = cost_old;  b = b_old;  xhat = xhat_old;
    c = cost(x{it}, u{it}, x_d, u_d, Q, R, Qf, num_steps);
    b = b_old;
    xhat = xhat_old;
    % x{it+1} = x{it};
    % u{it+1} = u{it};
%     s=s_old; Sinv=Sinv_old; v=v_old; Vinv=Vinv_old; r=r_old; R=R_old;
    s = s_old;
    Sinv = Sinv_old;
    v = v_old;
    Vinv = Vinv_old;
%   }else{
    else
%     //cout <<"\b AICO ACCEPT" <<endl;
%     damping /= 5.;
    Damping = Damping / 5;
%   }
    end
end

if break_when
    if it > 2 && abs(c - cost(x{it-1}, u{it-1}, x_d, u_d, Q, R, Qf, num_steps)) < thresh
        msg = ['Cost:', num2str(c)];
        disp(msg)
        num_iter = it
        break
    end
end

if pendulum
    msg = ['Cost:', num2str(c)];
    disp(msg)
    % tIter = cputime - tIter;
    % disp(tIter)
    % tIter = cputime;
    % msg = ['Theta Error:', num2str(sum(abs(x_d(1, :)-x{it+1}(1, :))))];
    % disp(msg);
    % msg = ['Final Theta Error:', num2str(x_d(1, end)-x{it+1}(1, end))];
    % disp(msg)
    % msg = ['Omega Error:', num2str(sum(abs(x_d(2, :)-x{it+1}(2, :))))];
    % disp(msg);
    % msg = ['Final Omega Error:', num2str(x_d(2, end)-x{it+1}(2, end))];
    % disp(msg)

    % msg = ['Control Error:', num2str(sum(abs(u_d(1, :)-u{it+1}(1, :))))];
    % disp(msg);
    % msg = ['Final Control Error:', num2str(u_d(1, end) - u{it+1}(1, end))];
    % disp(msg);

elseif cartpole
    msg = ['Position Error:', num2str(sum(abs(x_d(1, :)-x{it+1}(1, :))))];
    disp(msg);
    msg = ['Final Position Error:', num2str(x_d(1, end)-x{it+1}(1, end))];
    disp(msg)
    msg = ['Theta Error:', num2str(sum(abs(x_d(2, :)-x{it+1}(2, :))))];
    disp(msg);
    msg = ['Final Theta Error:', num2str(x_d(2, end)-x{it+1}(2, end))];
    disp(msg)
    msg = ['Velocity Error:', num2str(sum(abs(x_d(3, :)-x{it+1}(3, :))))];
    disp(msg);
    msg = ['Final Velocity Error:', num2str(x_d(3, end)-x{it+1}(3, end))];
    disp(msg)
    msg = ['Omega Error:', num2str(sum(abs(x_d(4, :)-x{it+1}(4, :))))];
    disp(msg);
    msg = ['Final Omega Error:', num2str(x_d(4, end)-x{it+1}(4, end))];
    disp(msg)

    msg = ['Control Error:', num2str(sum(abs(u_d(1, :)-u{it+1}(1, :))))];
    disp(msg);
    msg = ['Final Control Error:', num2str(u_d(1, end) - u{it+1}(1, end))];
    disp(msg)

elseif drones
    msg = ['x1 Error:', num2str(sum(abs(x_d(1, :)-x{it+1}(1, :))))];
    disp(msg);
    msg = ['Final x1 Error:', num2str(x_d(1, end)-x{it+1}(1, end))];
    disp(msg)
    msg = ['y1 Error:', num2str(sum(abs(x_d(2, :)-x{it+1}(2, :))))];
    disp(msg);
    msg = ['Final y1 Error:', num2str(x_d(2, end)-x{it+1}(2, end))];
    disp(msg)
    msg = ['x2 Error:', num2str(sum(abs(x_d(3, :)-x{it+1}(3, :))))];
    disp(msg);
    msg = ['Final x2 Error:', num2str(x_d(3, end)-x{it+1}(3, end))];
    disp(msg)
    msg = ['y2 Error:', num2str(sum(abs(x_d(4, :)-x{it+1}(4, :))))];
    disp(msg);
    msg = ['Final y2 Error:', num2str(x_d(4, end)-x{it+1}(4, end))];
    disp(msg)
    msg = ['x3 Error:', num2str(sum(abs(x_d(5, :)-x{it+1}(5, :))))];
    disp(msg);
    msg = ['Final x3 Error:', num2str(x_d(5, end)-x{it+1}(5, end))];
    disp(msg)
    msg = ['y3 Error:', num2str(sum(abs(x_d(6, :)-x{it+1}(6, :))))];
    disp(msg);
    msg = ['Final y3 Error:', num2str(x_d(6, end)-x{it+1}(6, end))];
    disp(msg)

    msg = ['Control x1 Error:', num2str(sum(abs(u_d(1, :)-u{it+1}(1, :))))];
    disp(msg);
    msg = ['Final Control x1 Error:', num2str(u_d(1, end) - u{it+1}(1, end))];
    disp(msg)
    msg = ['Control y1 Error:', num2str(sum(abs(u_d(2, :)-u{it+1}(2, :))))];
    disp(msg);
    msg = ['Final Control y1 Error:', num2str(u_d(2, end) - u{it+1}(2, end))];
    disp(msg)
    msg = ['Control x2 Error:', num2str(sum(abs(u_d(3, :)-u{it+1}(3, :))))];
    disp(msg);
    msg = ['Final Control x2 Error:', num2str(u_d(3, end) - u{it+1}(3, end))];
    disp(msg)
    msg = ['Control y2 Error:', num2str(sum(abs(u_d(4, :)-u{it+1}(4, :))))];
    disp(msg);
    msg = ['Final Control y2 Error:', num2str(u_d(4, end) - u{it+1}(4, end))];
    disp(msg)
    msg = ['Control x3 Error:', num2str(sum(abs(u_d(5, :)-u{it+1}(5, :))))];
    disp(msg);
    msg = ['Final Control x3 Error:', num2str(u_d(5, end) - u{it+1}(5, end))];
    disp(msg)
    msg = ['Control y3 Error:', num2str(sum(abs(u_d(6, :)-u{it+1}(6, :))))];
    disp(msg);
    msg = ['Final Control y3 Error:', num2str(u_d(6, end) - u{it+1}(6, end))];
    disp(msg)

end

end

tTotal = cputime - tStart;
disp(tTotal)
disp(tTotal/num_iter)

T = num_steps * dt;
time = linspace(0, T, num_steps);

if pendulum

figure(1)
hold on
for i = 1:num_iter+1
plot(x{i}(1, :), x{i}(2, :));
end
plot(x{num_iter+1}(1, :), x{num_iter+1}(2, :), 'LineWidth', 5);
plot(x_d(1, :), x_d(2, :), 'LineWidth', 2.5);
xlabel('theta (rad)');
ylabel('thetadot (rad/s)');

figure(2)
hold on
for i = 1:num_iter+1
plot(time, u{i}(1, :));
end
% plot(time, u{num_iter+1}(1, :));
plot(time, u_d(1, :));
xlabel('time (s)');
ylabel('control');

elseif cartpole

figure(1)
hold on
% for i = 1:num_iter+1
% plot(x{i}(1, :), x{i}(3, :));
% end
plot(x{num_iter+1}(1, :), x{num_iter+1}(3, :));
plot(x_d(1, :), x_d(3, :));
xlabel('position (m)');
ylabel('velocity (m/s)');

figure(2)
hold on
for i = 1:num_iter+1
plot(x{i}(2, :), x{i}(4, :));
end
% plot(x{num_iter+1}(2, :), x{num_iter+1}(4, :));
plot(x_d(2, :), x_d(4, :));
xlabel('theta (rad)');
ylabel('thetadot (rad/s)');

figure(3)
hold on
for i = 1:num_iter+1
plot(time, u{i}(1, :));
end
% plot(time, u{num_iter+1}(1, :));
plot(time, u_d(1, :));
xlabel('time (s)');
ylabel('control');

elseif drones

% figure(1)
% hold off
% % D1 = rectangle('Position', [x{num_iter+1}(1, 1), x{num_iter+1}(2, 1), 1, 1], 'Curvature', [1, 1]);
% % D2 = rectangle('Position', [x{num_iter+1}(3, 1), x{num_iter+1}(4, 1), 1, 1], 'Curvature', [1, 1]);
% % D3 = rectangle('Position', [x{num_iter+1}(5, 1), x{num_iter+1}(6, 1), 1, 1], 'Curvature', [1, 1]);
% % axis([-10, 10, -10, 10])
% % for i = 1:num_steps
% %     set(D1, 'Position', [x{num_iter+1}(1, i), x{num_iter+1}(2, i), 1, 1]);
% %     set(D2, 'Position', [x{num_iter+1}(3, i), x{num_iter+1}(4, i), 1, 1]);
% %     set(D3, 'Position', [x{num_iter+1}(5, i), x{num_iter+1}(6, i), 1, 1]);
% %     S12 = line([x{num_iter+1}(1, i), x{num_iter+1}(2, i)], [x{num_iter+1}(3, i), x{num_iter+1}(4, i)]);
% %     S23 = line([x{num_iter+1}(3, i), x{num_iter+1}(4, i)], [x{num_iter+1}(5, i), x{num_iter+1}(6, i)]);
% %     F(i) = getframe;
% % end
% D1 = rectangle('Position', [x_d(1, 1), x_d(2, 1), 1, 1], 'Curvature', [1, 1]);
% D2 = rectangle('Position', [x_d(3, 1), x_d(4, 1), 1, 1], 'Curvature', [1, 1]);
% D3 = rectangle('Position', [x_d(5, 1), x_d(6, 1), 1, 1], 'Curvature', [1, 1]);
% axis([-10, 10, -10, 10])
% for i = 1:num_steps
%     F(i) = getframe;
%     set(D1, 'Position', [x_d(1, i), x_d(2, i), 1, 1]);
%     set(D2, 'Position', [x_d(3, i), x_d(4, i), 1, 1]);
%     set(D3, 'Position', [x_d(5, i), x_d(6, i), 1, 1]);
%     S12 = line([x_d(1, i), x_d(2, i)], [x_d(3, i), x_d(4, i)]);
%     S23 = line([x_d(3, i), x_d(4, i)], [x_d(5, i), x_d(6, i)]);
% end
% F(end+1) = getframe;
% movie(F, 10)

figure(2)
hold on
plot(x_d(1, :), x_d(2, :))
plot(x_d(3, :), x_d(4, :))
plot(x_d(5, :), x_d(6, :))
plot(x{num_iter+1}(1, :), x{num_iter+1}(2, :))
plot(x{num_iter+1}(3, :), x{num_iter+1}(4, :))
plot(x{num_iter+1}(5, :), x{num_iter+1}(6, :))
legend('d1d', 'd2d', 'd3d', 'd1', 'd2', 'd3')

% figure(3)
% hold off
% D1 = rectangle('Position', [x{num_iter+1}(1, 1), x{num_iter+1}(2, 1), 1, 1], 'Curvature', [1, 1]);
% D2 = rectangle('Position', [x{num_iter+1}(3, 1), x{num_iter+1}(4, 1), 1, 1], 'Curvature', [1, 1]);
% D3 = rectangle('Position', [x{num_iter+1}(5, 1), x{num_iter+1}(6, 1), 1, 1], 'Curvature', [1, 1]);
% axis([-10, 10, -10, 10])
% for i = 1:num_steps
%     F(i) = getframe;
%     set(D1, 'Position', [x{num_iter+1}(1, i), x{num_iter+1}(2, i), 1, 1]);
%     set(D2, 'Position', [x{num_iter+1}(3, i), x{num_iter+1}(4, i), 1, 1]);
%     set(D3, 'Position', [x{num_iter+1}(5, i), x{num_iter+1}(6, i), 1, 1]);
%     S12 = line([x{num_iter+1}(1, i), x{num_iter+1}(2, i)], [x{num_iter+1}(3, i), x{num_iter+1}(4, i)]);
%     S23 = line([x{num_iter+1}(3, i), x{num_iter+1}(4, i)], [x{num_iter+1}(5, i), x{num_iter+1}(6, i)]);
% end
% F(end+1) = getframe;
% movie(F, 10)

end
