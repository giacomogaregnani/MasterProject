function [Bq, Fq, micMesh, micVp, micFemspace] = get_parameter_independent_forms(micMesh, micVp, Q)

Bq = cell(Q,1);
Fq = cell(Q,1);
micVpXt = micVp;
for q = 1:Q
    micVpXt.a = @(x,kk,ll) feval(micVp.a,x,kk,ll,q);
    [Bq{q}, Fq{q}, micMesh, micVp_tmp, micFemspace] = get_forms(micMesh, micVpXt);
end

a_tmp = micVp.a;
micVp = micVp_tmp;
micVp.a = a_tmp;