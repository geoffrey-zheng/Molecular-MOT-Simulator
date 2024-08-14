#0) Specify whether running on Macbook or Windows desktop
computer_type = 0 #0 for desktop, 1 for Macbook

#1)Go to directory and load external variables + functions
cd(@__DIR__);#moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well
include("AgVariables.jl") #change this to whatever molecule you care about
include("auxFunctionsCoinage.jl");#supplementary functions

#2) User choices with respect to saving output files
saveInRealUnits = 1;#if 1, save vel+accel in m/s, mm/ms^2.  If 0, save in normalized units (vel= v/(gam/k)), (force=1e-3*hbar*k*gam)
saveData=1;#if you want to save the data
saveDataFolderTag = "AgRedMOT"; #If you want to put anything additional in "savefoldername" to tag it, see variable folderString after lasers are declared.
addHeaders=1;

#3) Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
bGradReal = 10;# in units Gauss/cm.  if "Static", this becomes the static field, in Gauss
bGrad = (1 / kA * 1e2) * bGradReal;
numTrialsPerValueSet = 2;#number of trials per set of values (displacementsInMM,userSpeeds,longSpeeds)
velDirRelToR = 0;#-1 will choose random values for direction of v,r.  0 will force them to be same direction. 1 forces orthogonal.  2 forces opposite.
forceXY=1; #if '1', will force v, r to go along (x+y)/sqrt(2).  Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = 0
if velDirRelToR==-1
    runType = string("Random");#goes into folder name.
elseif velDirRelToR==0
    runType = string("SameDir");#goes into folder name.
elseif velDirRelToR==1
    runType = string("OrthoDir");#goes into folder name.
else
    runType = string("OppDir");#goes into folder name.
end
vRound = 0.02;#nearest normalized unit to round all speeds to.  Simulation will have periodicity of 2*pi/vRound, so if you round finely, will take a while.  Usually choose 0.02

#4A) parameters for quick test of restoring force
#=
displacementsInMM = [1,2,3,5,7,9,11,14,17];
userSpeeds = [-8,-6.5,-5,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,5,6.5,8];
=#
displacementsInMM = [0.5,1.5,3.0,4.5,6.0,7.5,9.0,10.5,12,13.5,15];
userSpeeds = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-.8,-.6,-.4,-.2,-.1,-.05,.05,.1,.2,.4,.6,.8,1,1.5,2,2.5,3,3.5,4];
#displacementsInMM = [0.2,0.5,1,2,3,4,5];
#displacementsInMM = [0.1,0.2,0.3,0.5,0.7,1,1.5,2];
#userSpeeds = [-3,-1,-.3,.3,1,3]
#blueMOT
#userSpeeds = [-5,-4,-3,-2.5,-2,-1.5,-1,-.5,-.3,-0.2,-0.1,-0.05,-.03,-0.02,-0.01,.01,.02,.03,0.05,0.1,0.2,.3,.5,1,1.5,2,2.5,3,4,5]
#userSpeeds = [-1,-.5,-.3,-.1,-.05,-.02,-.01,-0.005,-0.003,0.003,0.005,.01,.02,.05,.1,.3,.5,1];
#displacementsInMM = [.5,1,1.5];
#userSpeeds = [-3,-2.5,-2,-1.6,-1.4,-1.2,-1,-.8,-.6,-.5,-.4,-.3,-.25,-.2,-.15,-.1,-.07,-.05,-.03,.03,.05,.07,.1,.15,.2,.25,.3,.4,.5,.6,.8,1,1.2,1.4,1.6,2,2.5,3];

#laser parameters
#=
s0 = [1.0,1.0,1.0].*1.0;
laserEnergy = [-15.8,7.2,-14.8];
polSign = [-1,-1,1];
whichJ = ["1/2","1/2","1/2"];
whichFGround = ["FHigh","FLow","FHigh"]
waists = [5,5,5] .* 1e-3 .* kA;#in normalized units
=#
s0 = [.8,.1].*1.0;
laserEnergy = [-3.0,0.];
polSign = [-1,-1];
whichJ = ["3/2","3/2"];
whichFGround = ["FHigh","FLow"]
waists = [5,5] .* 1e-3 .* kA;#in normalized units, converted from mm

(couplingMatrices,bCouplingMatrices,hyperfineMatrix,stateEnergyMatrix,laserMasks,numZeemanStatesGround,numZeemanStatesExcited) = createCouplingTermsandLaserMasks(whichJ,whichFGround)

lasers = Lasers(s0,laserEnergy,polSign,whichJ,whichFGround,waists,laserMasks);#define lasers structure, see auxFunctions

numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited;

#everything below here you is part of the initizliation block that you probably won't want to change

rInit = [0.0, 0.0, 0];#placehold not used
vInit = [3.0, 0.0, 0.1];#placehold not used
rpInit = zeros(length(s0),6);
#note: p will include a lot of pre-allocated stuff.  This is basically all of the stuff 'passed to' the obe solver, in addition to the initial condition of the density matrix defined below
#in retrospect p is not the best choice for the variable name but it's the julia house style...maybe replace later. (actually you can't.  Julia forces ODEProblem to have a variable 'p')
pPreInitialized = preInitializer(length(s0),numZeemanStatesGround,numZeemanStatesTotal)


p = [rInit, vInit, rpInit, stateEnergyMatrix, lasers, bGrad * normalizedBohrMag,
    couplingMatrices[1], couplingMatrices[2], couplingMatrices[3], bCouplingMatrices[1], bCouplingMatrices[2], bCouplingMatrices[3],hyperfineMatrix];
append!(p,pPreInitialized);

#density matrix Initialization
pStart = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
pStart[1,1] = 1/2;
pStart[2,2] = 1/2;

bichrom=0;
oneHalf=0;
repump=0;
if "1/2" in whichJ
    if "3/2" in whichJ
        bichrom=1;
    else
        oneHalf = 1;
    end
end
if "FLow" in whichFGround
    repump=1;
end


#initialize some 'masks' that zero out subset of population values...helpful for quick calculation of populations in various ground states
nFLow = trunc(Int,2*(1/2+nucSpin-1)+1);
nFHigh = trunc(Int,2*(1/2+nucSpin)+1);
nFeLowOneHalf = trunc(Int,2*(1/2+nucSpin-1)+1);
nFeHighOneHalf = trunc(Int,2*(1/2+nucSpin)+1);
nFeLowThreeHalf = trunc(Int,2*(3/2+nucSpin-1)+1);
nFeHighThreeHalf = trunc(Int,2*(3/2+nucSpin)+1);

maskFHigh = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFLow = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFeOneHalfHigh = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFeOneHalfLow = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFeThreeHalfHigh = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFeThreeHalfLow = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
pHighVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pLowVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pHighOneHalfVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pLowOneHalfVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pHighThreeHalfVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pLowThreeHalfVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);

if repump==1
    if bichrom==1  
        headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFHigh" "PFLow" "PFeOneHalfHigh"  "PFeOneHalfLow" "PFeThreeHalfHigh" "PFeThreeHalfLow"];

        maskFLow[1:nFLow, 1:nFLow] .= ones(nFLow, nFLow);
        maskFHigh[(nFLow+1):(nFLow+nFHigh),(nFLow+1):(nFLow+nFHigh)] .= ones(nFHigh, nFHigh);

        maskFeOneHalfLow[(nFLow+nFHigh+1):(nFLow+nFHigh+nFeLowOneHalf),(nFLow+nFHigh+1):(nFLow+nFHigh+nFeLowOneHalf)] .= ones(nFeLowOneHalf, nFeLowOneHalf);
        maskFeOneHalfHigh[(nFLow+nFHigh+nFeLowOneHalf+1):(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf),(nFLow+nFHigh+nFeLowOneHalf+1):(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf)] .= ones(nFeHighOneHalf, nFeHighOneHalf);

        maskFeThreeHalfLow[(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf+1):(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf+nFeLowThreeHalf),(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf+1):(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf+nFeLowThreeHalf)] .= ones(nFeLowThreeHalf, nFeLowThreeHalf);
        maskFeThreeHalfHigh[(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf+nFeLowThreeHalf+1):numZeemanStatesTotal,(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf+nFeLowThreeHalf+1):numZeemanStatesTotal] .= ones(nFeHighThreeHalf, nFeHighThreeHalf);
    elseif oneHalf ==1
        headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFHigh" "PFLow" "PFeOneHalfHigh"  "PFeOneHalfLow"];

        maskFLow[1:nFLow, 1:nFLow] .= ones(nFLow, nFLow);
        maskFHigh[(nFLow+1):(nFLow+nFHigh),(nFLow+1):(nFLow+nFHigh)] .= ones(nFHigh, nFHigh);

        maskFeOneHalfLow[(nFLow+nFHigh+1):(nFLow+nFHigh+nFeLowOneHalf),(nFLow+nFHigh+1):(nFLow+nFHigh+nFeLowOneHalf)] .= ones(nFeLowOneHalf, nFeLowOneHalf);
        maskFeOneHalfHigh[(nFLow+nFHigh+nFeLowOneHalf+1):(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf),(nFLow+nFHigh+nFeLowOneHalf+1):(nFLow+nFHigh+nFeLowOneHalf+nFeHighOneHalf)] .= ones(nFeHighOneHalf, nFeHighOneHalf);

    else
        headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFHigh" "PFLow" "PFeThreeHalfHigh" "PFeThreeHalfLow"];
        maskFLow[1:nFLow, 1:nFLow] .= ones(nFLow, nFLow);
        maskFHigh[(nFLow+1):(nFLow+nFHigh),(nFLow+1):(nFLow+nFHigh)] .= ones(nFHigh, nFHigh);

        maskFeThreeHalfLow[(nFLow+nFHigh+1):(nFLow+nFHigh+nFeLowThreeHalf),(nFLow+nFHigh+1):(nFLow+nFHigh+nFeLowThreeHalf)] .= ones(nFeLowThreeHalf, nFeLowThreeHalf);
        maskFeThreeHalfHigh[(nFLow+nFHigh+nFeLowThreeHalf+1):numZeemanStatesTotal,(nFLow+nFHigh+nFeLowThreeHalf+1):numZeemanStatesTotal] .= ones(nFeHighThreeHalf, nFeHighThreeHalf);
    end

else
    if bichrom==1  
        headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFHigh" "PFeOneHalfHigh" "PFeThreeHalfHigh"];

        maskFHigh[(1):(nFHigh),(1):(nFHigh)] .= ones(nFHigh, nFHigh);

        maskFeOneHalfHigh[(nFHigh+1):(nFHigh+nFeHighOneHalf),(nFHigh+1):(nFHigh+nFeHighOneHalf)] .= ones(nFeHighOneHalf, nFeHighOneHalf);

        maskFeThreeHalfHigh[(nFHigh+nFeHighOneHalf+1):numZeemanStatesTotal,(nFHigh+nFeHighOneHalf+1):numZeemanStatesTotal] .= ones(nFeHighThreeHalf, nFeHighThreeHalf);
    elseif oneHalf ==1
        headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFHigh" "PFeOneHalfHigh"];
        maskFHigh[(1):(nFHigh),(1):(nFHigh)] .= ones(nFHigh, nFHigh);

        maskFeOneHalfHigh[(nFHigh+1):(nFHigh+nFeHighOneHalf),(nFHigh+1):(nFHigh+nFeHighOneHalf)] .= ones(nFeHighOneHalf, nFeHighOneHalf);

    else
        headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFHigh" "PFeThreeHalfHigh"];

        maskFHigh[(1):(nFHigh),(1):(nFHigh)] .= ones(nFHigh, nFHigh);

        maskFeThreeHalfHigh[(nFHigh+1):numZeemanStatesTotal,(nFHigh+1):numZeemanStatesTotal] .= ones(nFeHighThreeHalf, nFeHighThreeHalf);
    end
end

#save data folder path name
if saveData==1
    if computer_type == 0
        folderString = string(@__DIR__,"\\saveData\\",saveDataFolderTag, "BGradGPerCM",bGradReal,"Date",Dates.format(now(),"yyyymmdd_HHMM"))
        mkpath(folderString)
    else
        folderString = string(@__DIR__,"/saveData/",saveDataFolderTag, "BGradGPerCM",bGradReal,"Date",Dates.format(now(),"yyyymmdd_HHMM"))
        mkpath(folderString)
    end
end

#below still needs to be edited (TO DO)


forceVsTime = Array{Array{ComplexF64,2},1}(undef, numTrialsPerValueSet * 2);
forceVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot v/|v|
forceVsPos = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot r/|r|

#8) Iterate over user choices for displacements and speeds

for l = 1:length(displacementsInMM)
    currDisp = displacementsInMM[l];
    for j = 1:length(userSpeeds)
        currSpeed = userSpeeds[j];
        if abs(currSpeed)<0.01
            vRound = 0.001;
        elseif abs(currSpeed)<0.03
            vRound = 0.002;
        elseif abs(currSpeed)<0.1
            vRound = 0.01;
        elseif abs(currSpeed)<0.5
            vRound = 0.02;
        else
            vRound = 0.05;
        end
        #8A) Set up and solve OBEs
        (randRxs,randRys,randRzs,randVxs,randVys,randVzs,rp) = generateRandPosAndVel(numTrialsPerValueSet,velDirRelToR,currDisp,currSpeed,vRound,forceXY,length(s0));
        tForSteadyState = maximum([10 / currSpeed, 1800]);#obtained by trial and error.  Could potentially be handled more rigrorously (solve ode in steps of 'period length' until solution 'converges')
        periodLength = 2 * pi / vRound;
        saveTimes = tForSteadyState:0.1:(tForSteadyState+periodLength)#times to record obe solution for force integration
        for i = 1:(numTrialsPerValueSet*2)
            forceVsTime[i] = zeros(length(saveTimes), 3)
        end
        prob = ODEProblem(densityMatrixChangeTerms!, pStart, (0.0, tForSteadyState + periodLength), p)#set up OBE problem to solve

        function prob_func(prob, i, repeat)#change position and velocity of funtion each 'ensemble' sample based on the randomly chosen values for pos/vel vector (magnitude fixed)
            prob.p[1][1] = randRxs[i]
            prob.p[1][2] = randRys[i]
            prob.p[1][3] = randRzs[i]
            prob.p[2][1] = randVxs[i]
            prob.p[2][2] = randVys[i]
            prob.p[2][3] = randVzs[i]
            #=
            for laserVar = 1:length(s0)
                prob.p[3][laserVar,1] = rp[laserVar,1] .* useRandPhase;
                prob.p[3][laserVar,2] = rp[laserVar,2] .* useRandPhase;
                prob.p[3][laserVar,3] = rp[laserVar,3] .* useRandPhase;
                prob.p[3][laserVar,4] = rp[laserVar,4] .* useRandPhase;
                prob.p[3][laserVar,5] = rp[laserVar,5] .* useRandPhase;
                prob.p[3][laserVar,6] = rp[laserVar,6] .* useRandPhase;
            end
            =#
            remake(prob)
        end

        #these two lines here actually handle the parallized runs of the ode solver
        ens_prob = EnsembleProblem(prob, prob_func=prob_func)#solve obe problem for various initial conditions re-set by 'prob_func' each iteration
        @time sol = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=numTrialsPerValueSet * 2, saveat=saveTimes)#parallelized OBE solver, runs on amount of threads made available by CPU (Threads.nthreads())
        
        #8B) calculate forces (f\dot r/|r|, etc.) for each random R, V trial..
        @time for i = 1:(numTrialsPerValueSet*2)
            currSol = sol[i]
            makeForceVsTime!(forceVsTime[i], currSol.t, currSol.u, lasers,
            couplingMatrices, [randRxs[i], randRys[i], randRzs[i]], [randVxs[i], randVys[i], randVzs[i]],rp)

            forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
            randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
            randVzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) / (currSol.t[end] - currSol.t[1])
            forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
            randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
            randRzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) / (currSol.t[end] - currSol.t[1])

            pHighVsSpeed[j, i] = mean(real(tr.([maskFHigh .* v for v in currSol.u])))
            pLowVsSpeed[j, i] = mean(real(tr.([maskFLow .* v for v in currSol.u])))
            pHighOneHalfVsSpeed[j, i] = mean(real(tr.([maskFeOneHalfHigh .* v for v in currSol.u])))
            pLowOneHalfVsSpeed[j, i] = mean(real(tr.([maskFeOneHalfLow .* v for v in currSol.u])))
            pHighThreeHalfVsSpeed[j, i] = mean(real(tr.([maskFeThreeHalfHigh .* v for v in currSol.u])))
            pLowThreeHalfVsSpeed[j, i] = mean(real(tr.([maskFeThreeHalfLow .* v for v in currSol.u])))
        end#for all trials
    end#for speeds

    #8C) for given set of speeds, for current choices of longSpeed and displacement, average a\dot v, a\dot r, populations,etc. over runs
    forceVsSpeedAvg = mean(forceVsSpeed, dims=2)
    forceVsSpeedAvg = dropdims(forceVsSpeedAvg, dims=(2))#converts to vector
    forceVsSpeedUnc = std(forceVsSpeed, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
    forceVsSpeedUnc = dropdims(forceVsSpeedUnc, dims=(2))

    forceVsPosAvg = mean(forceVsPos, dims=2)
    forceVsPosAvg = dropdims(forceVsPosAvg, dims=(2))
    forceVsPosUnc = std(forceVsPos, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
    forceVsPosUnc = dropdims(forceVsPosUnc, dims=(2))

    pHighVsSpeedAvg = mean(pHighVsSpeed, dims=2)
    pHighVsSpeedAvg = dropdims(pHighVsSpeedAvg, dims=(2))
    pLowVsSpeedAvg = mean(pLowVsSpeed, dims=2)
    pLowVsSpeedAvg = dropdims(pLowVsSpeedAvg, dims=(2))
    pHighOneHalfVsSpeedAvg = mean(pHighOneHalfVsSpeed, dims=2)
    pHighOneHalfVsSpeedAvg = dropdims(pHighOneHalfVsSpeedAvg, dims=(2))
    pLowOneHalfVsSpeedAvg = mean(pLowOneHalfVsSpeed, dims=2)
    pLowOneHalfVsSpeedAvg = dropdims(pLowOneHalfVsSpeedAvg, dims=(2))
    pHighThreeHalfVsSpeedAvg = mean(pHighThreeHalfVsSpeed, dims=2)
    pHighThreeHalfVsSpeedAvg = dropdims(pHighThreeHalfVsSpeedAvg, dims=(2))
    pLowThreeHalfVsSpeedAvg = mean(pLowThreeHalfVsSpeed, dims=2)
    pLowThreeHalfVsSpeedAvg = dropdims(pLowThreeHalfVsSpeedAvg, dims=(2))

    #8D) convert to real units if applicable and save data
    (forceVsSpeedAvgSaveVals,forceVsSpeedUncSaveVals,forceVsPosAvgSaveVals,forceVsPosUncSaveVals) = (forceVsSpeedAvg,forceVsSpeedUnc,forceVsPosAvg,forceVsPosUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
    userSpeedsSaveVals = userSpeeds.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));
    if saveData==1
        if computer_type == 0
            open(string(folderString,"/forceVsSpeedDisplacement",displacementsInMM[l],"MM",runType,".dat"),"a") do io
                if addHeaders==1
                    if repump == 1
                        if bichrom==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg)])
                        elseif oneHalf==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg)])
                        else
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg)])
                        end
                    else
                        if bichrom==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg)])
                        elseif oneHalf==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg)])
                        else
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighThreeHalfVsSpeedAvg)])
                        end
                    end
                        

                else #if you've already added headers/don't want them, just append the current forceVsSpeed to the relevant file (so, if you have different longSpeeds, they'll all show up in same file since file is distinguished by displacement)
                    if repump == 1
                        if bichrom==1
                            writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg))
                        elseif oneHalf==1
                            writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg))
                        else
                            writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg))
                        end
                    else
                        if bichrom==1
                            writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg))
                        elseif oneHalf==1
                            writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg))
                        else
                            writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighThreeHalfVsSpeedAvg))
                        end
                    end
                end
            end
        else
            open(string(folderString,"/forceVsSpeedDisplacement",displacementsInMM[l],"MM",runType,".dat"),"a") do io
                if addHeaders==1
                    if repump == 1
                        if bichrom==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg)])
                        elseif oneHalf==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg)])
                        else
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg)])
                        end
                    else
                        if bichrom==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg)])
                        elseif oneHalf==1
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg)])
                        else
                            writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighThreeHalfVsSpeedAvg)])
                        end
                    end
                        

                else #if you've already added headers/don't want them, just append the current forceVsSpeed to the relevant file (so, if you have different longSpeeds, they'll all show up in same file since file is distinguished by displacement)
                    if repump == 1
                        if bichrom==1
                            writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg))
                        elseif oneHalf==1
                            writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighOneHalfVsSpeedAvg, pLowOneHalfVsSpeedAvg))
                        else
                            writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pLowVsSpeedAvg, pHighThreeHalfVsSpeedAvg,pLowThreeHalfVsSpeedAvg))
                        end
                    else
                        if bichrom==1
                            writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg, pHighThreeHalfVsSpeedAvg))
                        elseif oneHalf==1
                            writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighOneHalfVsSpeedAvg))
                        else
                            writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pHighVsSpeedAvg, pHighThreeHalfVsSpeedAvg))
                        end
                    end
                end
            end
        end
    end
end#for displacements

#save laser variables used to file
laserVarHeaders = ["s0" "energy" "polSign" "whichJ" "whichFGround"]
if saveData ==1
    if computer_type == 0
        open(string(folderString,"\\laserVariables.dat"),"w") do io
            writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign,whichJ,whichFGround)]);
        end
    else
        open(string(folderString,"/laserVariables.dat"),"w") do io
            writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign,whichJ,whichFGround)]);
        end
    end
end

