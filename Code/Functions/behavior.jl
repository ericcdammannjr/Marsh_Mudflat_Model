function behavior()

    function condition(out, u, t, integrator)

        out[1] = u[3]-(0.7167*r-0.483)  # Drown
        out[2] = 0.99*L-u[1]         # Erode
        out[3] = u[1]-0.01*L         # Fill

    end

    function affect!(integrator, idx)

        if idx == 1

            terminate!(integrator)

        elseif idx == 2

            terminate!(integrator)

        elseif idx == 3

            terminate!(integrator)

        end

    end

    cb = VectorContinuousCallback(condition,affect!,3)

    return cb

end