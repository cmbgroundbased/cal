module CAL
    using CxxWrap
    @wrapmodule "build/lib/libAtmSimJulia"

    """
        Copy the python docstring

    """
    function atm_get_abs_coef(altitude::Float64, temperature::Float64, pressure::Float64, pwv::Float64, freq::Float64)

        return CAL.atm_get_absorption_coefficient(altitude, temperature, pressure, pwv, freq)

    end

    function atm_get_abs_coef_vec(altitude::Float64, temperature::Float64, pressure::Float64, pwv::Float64, freq_min::Float64, freq_max::Float64, nfreq::Int64)

        Buff = Array{Float64, 1}(undef, nfreq)
        CAL.atm_get_absorption_coefficient_vec(altitude, temperature, pressure, pwv, freq_min, freq_max, nfreq, Buff)

        return Buff
    end

    function atm_get_atm_load(altitude::Float64, temperature::Float64, pressure::Float64, pwv::Float64, freq::Float64)

        return CAL.atm_get_atmospheric_loading(altitude, temperature, pressure, pwv, freq)

    end

    function atm_get_atm_load_vec(altitude::Float64, temperature::Float64, pressure::Float64, pwv::Float64, freq_min::Float64, freq_max::Float64, nfreq::Int64)

        Buff = Array{Float64, 1}(undef, nfreq)
        CAL.atm_get_atmospheric_loading_vec(altitude, temperature, pressure, pwv, freq_min, freq_max, nfreq, Buff)
        return Buff

    end


    function __init__()
        @initcxx
    end

end


