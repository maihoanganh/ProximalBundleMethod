using ProximalBundleMethod

"""
min f(x) x in R^n

f is nonsmooth function
g is a subgradient of f
"""


n=100

function evaluate_f(x::Vector{Float64})

  f = 0.0;
  g=zeros(Float64,n)

  for i=1:n-1
    g[i+1] = 0.0;
    a = -x[i]-x[i+1];
    b = -x[i]-x[i+1]+(x[i]*x[i]+x[i+1]*x[i+1]-1.0);

    if (a >= b)
      f = f+a;
      g[i] = g[i]-1.0;
      g[i+1] = -1.0;
    else
      f = f+b;  
      g[i] = g[i]-1.0+2.0*x[i];
      g[i+1] = -1.0+2*x[i+1];
    end
  end
    return f, g
end

tol=1e-3

bundle = ProximalBundleMethod.Model{ProximalBundleMethod.ProximalMethod}(n, evaluate_f,tol)

ProximalBundleMethod.runb(bundle)

optsol=ProximalBundleMethod.getsolution(bundle)
optval=ProximalBundleMethod.getobjectivevalue(bundle)
