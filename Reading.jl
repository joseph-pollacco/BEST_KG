
module reading
    export  VANGENUCHTEN, INFILTRATION, BEST, HYDRAULIC

    #= =============== READING INFILTRATION DATA =============== =#
    function INFILTRATION(Path_best)
        Data = Array(readdlm(Path_best, ',', header=true, skipstart=0)[1]) # skip first line with header
        ID_1 = Data[:,1]  # Reading column 1
        T_1 = Data[:,2] # Reading column 2
        Infiltration_1 =  Data[:,3] # Reading column 3
        Data_N = length(ID_1) # Total number of data
        Sample_N = Int(maximum(ID_1[:])) # Total number of soil samples since the ID is writen in order of the lowest to highest

        #= Reads the infiltration data in column and the data transformed to Array[iS, :] based on the id of the soil sample =#
        ID = Array{Float64}(Sample_N , 10000)
        Time = Array{Float64}(Sample_N , 10000)
        Infiltration = Array{Float64}(Sample_N , 10000)
        Time_Ni = Array{Int}(Sample_N) # Number of infiltration points for each samples
        iS = 1
        i_Infiltration = 0
        ID_Prev = 1. # value of ID(i-1)
        for i in 1:Data_N
            # IF we change the soil sample
             if ID_1[i] > ID_Prev
                iS +=1
                i_Infiltration = 1 #Starting from 1
            else
                iS +=0
                i_Infiltration += 1
            end
            ID_Prev = ID_1[i] #Previous ID number
            Time_Ni[iS] = i_Infiltration # Number of measurements in each soil sample,
            Time[iS, i_Infiltration] = T_1[i]
            Infiltration[iS, i_Infiltration] = Infiltration_1[i]
        end
        return Time, Infiltration, Time_Ni, Sample_N
    end


        #= =============== READING VANGENUCHTEN HYDRAULIC DATA =============== =#
     #= If available reading the van genuchten hydraulic params  =#
     function HYDRAULIC(Path)
        Data = Array(readdlm(Path, ',', header=true, skipstart=0)[1]) # skip first line with header\
        ID = Data[:,1]
        Se_Ini= Data[:,2]
        θs= Data[:,3]
        θr= Data[:,4]
        σ= Data[:,5]
        Hkg= Data[:,6]
        Ks= Data[:,7]
        Hvg= Data[:,8]
        N= Data[:,9]
        Km= Data[:,10]
        RingRadius= Data[:,11]
        SampleTrue= Int64.(Data[:,12])
        Name = string.(Data[:,13])
        return Se_Ini, θs, θr, σ, Hkg, Ks, Hvg, N, Km, RingRadius, SampleTrue, Name
    end

    
    #= =============== READING VANGENUCHTEN HYDRAULIC DATA =============== =#
     #= If available reading the van genuchten hydraulic params  =#
    function VANGENUCHTEN(Path_vg)
        Data = Array(readdlm(Path_vg, ',', header=true, skipstart=0)[1]) # skip first line with header\
        ID = Data[:,1]
        θr = Data[:,2]
        θs =  Data[:,3]
        Hvg = Data[:,4]
        n = Data[:,5]
        Km = Data[:,6]
        Ks = Data[:,7]
        return θs, θr, Hvg, n, Km, Ks
    end


 #= =============== READING BEST DATA WHICH DO NOT VARY WITH TIME =============== =#
    function BEST(Path_BEST)
        Data = Array(readdlm(Path_BEST, ',', header=true, skipstart=0)[1]) # skip first line with header\
        ID = Data[:,1]
        H_Ini = Data[:,2]
        RingRadius = Data[:,3]
        return H_Ini, RingRadius
    end
end
