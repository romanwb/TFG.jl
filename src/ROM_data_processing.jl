using MAT
#Lee el archivo original .mat y obtiene la informacion para el CB
function read_rom_input(path::String)
    matfile = matopen(path)
    Mee = read(matfile, "Mee")
    Mec = read(matfile, "Mec")
    Mcc = read(matfile, "Mcc")
    Kee = read(matfile, "Kee")
    Kec = read(matfile, "Kec")
    Kcc = read(matfile, "Kcc")
    F_vector = read(matfile, "F_vector")
    Rc = read(matfile, "Rc")
    close(matfile)

    # Ajusta los índices según el tamaño de fe y fc en tu modelo
    fe = F_vector[1:end-(size(Kcc,1))]
    fc = F_vector[end-(size(Kcc,1))+1:end]

    return (
        Mee=Mee, Mec=Mec, Mcc=Mcc,
        Kee=Kee, Kec=Kec, Kcc=Kcc,
        fe=fe, fc=fc, Rc=Rc
    )
end

# Con input el output de la funcion anterior, aplica el CB; como output: NamedTuple
function reduce_rom(data::NamedTuple; n_modes::Int)
    M, Max, Mxx, K, Kxx, ω², fₐ, fₓ, Φ, Ψ, T_CB = craig_bampton_aplication(
        data.Mee, data.Mec, data.Mcc,
        data.Kee, data.Kec, data.Kcc,
        data.fe, data.fc; n_modes=n_modes
    )

    return (
        M=M, K=K, Max=Max, Mxx=Mxx, Kxx=Kxx,
        omega2=ω², fa=fₐ, fx=fₓ, autovecs=Φ, psi=Ψ, T_CB=T_CB,
        n_modes=n_modes, Rc=data.Rc
    )
end

#Con input path donde escribirá el nuevo .mat con las matrices del CB, y tambien con input el output NamedTuple de la funcion anterior
function save_rom_data(path::String, reduced_data::NamedTuple)
    file = matopen(path, "w")

    for (key, value) in pairs(reduced_data)
        write(file, string(key), value)
    end

    close(file)
end

#Carga en el script los datos del nuevo .mat generado
function load_reduced_rom(path::String)
    file = matopen(path)
    data = (
        M       = real(read(file, "M")),
        K       = real(read(file, "K")),
        Max     = real(read(file, "Max")),
        Mxx     = real(read(file, "Mxx")),
        Kxx     = real(read(file, "Kxx")),
        ωₐ²     = real(read(file, "omega2")),
        fₐ      = real(read(file, "fa")),
        fₓ      = real(read(file, "fx")),
        Φ       = real(read(file, "autovecs")),
        Ψ       = real(read(file, "psi")),
        T_CB    = real(read(file, "T_CB")),
        n_modes = read(file, "n_modes"),
        Rₓ      = read(file, "Rc")
    )
    close(file)
    return data
end
