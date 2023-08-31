module MyLib
    #integral using trapezoid formula
    function IntTrap(x,y)
        N=length(x)
        ans=0.0
        ans+=2*sum(y)-y[N]-y[1]
        ans*=(x[N]-x[1])/(2*(N-1))
        return ans
    end

	# y[1]=y[i-2], y[2]=y[i-1], y[3]=y[i], y[4]=y[i+1], y[5]=y[i+2]
	function diff1st5pt(h,y)
		return (-y[5]/12 + y[4]*2/3 - y[2]*2/3 + y[1]/12)/h
	end

    function diff1st5pt(h,y1,y2,y3,y4,y5)
		return (-y5/12 + y4*2/3 - y2*2/3 + y1/12)/h
	end

	# only used for the rmesh=0:h:...
	# even Pality: P=1, odd Pality: P=-1
    function diff1st5pt(h::Float64,y::Vector{Float64},P::Int)
        N=length(y)
        dy=zeros(Float64,N)
        #dy[1]=(-y[3]/12 + y[2]*2/3 - P*y[2]*2/3 + P*y[3]/12)/h #rmesh=0:h:...
		dy[1]=(-y[3]/12 + y[2]*2/3 - P*y[1]*2/3 + P*y[2]/12)/h #rmesh=0.5*h:h:...
		#dy[2]=(-y[4]/12 + y[3]*2/3 - y[1]*2/3 + P*y[2]/12)/h #rmesh=0:h:...
        dy[2]=(-y[4]/12 + y[3]*2/3 - y[1]*2/3 + P*y[1]/12)/h #rmesh=0.5*h:h:...
		@views for i in 3:N-2
			dy[i]=diff1st5pt(h,y[i-2:i+2])
		end
		dy[N-1]=(y[N]-y[N-2])/(2*h)
		dy[N]=(-y[N-2]+4*y[N-1]-3*y[N])/(2*(-h))

        return dy
    end

	function diff2nd5pt(h,y)
		return (-y[5]/12 + y[4]*4/3 - y[3]*5/2 + y[2]*4/3 - y[1]/12)/h^2
	end

    function diff2nd5pt(h,y1,y2,y3,y4,y5)
		return (-y5/12 + y4*4/3 - y3*5/2 + y2*4/3 - y1/12)/h^2
	end

    function diff2nd5pt(h::Float64,y::Vector{Float64},P::Int)
        N=length(y)
        ddy=zeros(Float64,N)
		#ddy[1]=(-y[3]/12 + y[2]*4/3 - y[1]*5/2 + P*y[2]*4/3 - P*y[3]/12)/h^2 #rmesh=0:h:...
		ddy[1]=(-y[3]/12 + y[2]*4/3 - y[1]*5/2 + P*y[1]*4/3 - P*y[2]/12)/h^2 #rmesh=0.5*h:h:...
		#ddy[2]=(-y[4]/12 + y[3]*4/3 - y[2]*5/2 + y[1]*4/3 - P*y[2]/12)/h^2 #rmesh=0:h:...
		ddy[2]=(-y[4]/12 + y[3]*4/3 - y[2]*5/2 + y[1]*4/3 - P*y[1]/12)/h^2 #rmesh=0.5*h:h:...
		@views for i in 3:N-2
			ddy[i]=diff2nd5pt(h,y[i-2:i+2])
		end
		ddy[N-1]=(y[N]-2*y[N-1]+y[N-2])/(h^2)
		ddy[N]=(2*y[N]-5*y[N-1]+4*y[N-2]-y[N-3])/(h^2)

        return ddy
    end

    # ref. 計算物理学
    # solve Poisson eq. U''=-4πrρ(r)
    # qmax = \int d^3r ρ(r)
    function SolvePoissonEq(ρ::Vector{Float64},rmesh,qmax)
        N=length(rmesh)
        U=zeros(Float64,N)
        h=rmesh[2]-rmesh[1]

        #mesh for h:h:...
        #U[1]=rmesh[1]
        #U[2]=2*U[1]-h^2*4*π*rmesh[1]*ρ[1]

		#mesh for 0.5:h:...
		U[1]=rmesh[1]
		U[2]=3*U[1]-h^2*4*π*rmesh[1]*ρ[1]

        for i in 2:N-1
            U[i+1]=2*U[i]-U[i-1]-h^2*4*π*rmesh[i]*ρ[i]
        end

        #adjust the value by adding the solution of U''=0 (U=αr)
        α=(U[N]-qmax)/rmesh[N]
        @. U[:]-=α*rmesh[:]

        return U
    end

    #using isnan() to check the convergence
    function MyBisect(GivenLow, GivenUp, F::Function, args;rtol=1e-4)
        @assert GivenLow<GivenUp

        low=GivenLow
        up=GivenUp
        Flow=F(low,args...)
        Fup=F(up,args...)

        while abs((low-up)/(low+up)) > rtol
            if Flow*Fup>0
                return NaN
            else
                mean=(low+up)/2
                Fmean=F(mean,args...)
                if Fmean*Flow<0
                    up=mean
                else
                    low=mean
                end
            end
        end

        return (up+low)/2
    end

	const golden_ratio = (3-sqrt(5))/2
	# assuming F is heavy to calculate.
	# F should be include one argument.
	function GoldenSecSearch1Dim(GivenLow,GivenUp,F::Function;rtol=1e-4,maxiter=100)
		x_lower=GivenLow
		x_upper=GivenUp
		new_minimizer = x_lower + golden_ratio*(x_upper - x_lower)
		new_minimum = F(new_minimizer)
		count=0
		converged=false

		while count!=maxiter
			count+=1
			if (x_upper-x_lower)/(GivenUp-GivenLow)<rtol
				converged=true
				break
			end

			# lower minimizer            upper case
			if x_upper - new_minimizer > new_minimizer - x_lower
				new_x = new_minimizer + golden_ratio*(x_upper - new_minimizer)
				new_f = F(new_x)
				if new_f < new_minimum
					x_lower = new_minimizer
					new_minimizer = new_x
					new_minimum = new_f
				else
					x_upper = new_x
				end
			else # lower       minimizer upper case
				new_x = new_minimizer - golden_ratio*(new_minimizer - x_lower)
				new_f = F(new_x)
				if new_f < new_minimum
					x_upper = new_minimizer
					new_minimizer = new_x
					new_minimum = new_f
				else
					x_lower = new_x
				end
			end
		end

		if converged==false
			return NaN,NaN
		else
			return new_minimizer, new_minimum
		end
	end

end # module MyLib