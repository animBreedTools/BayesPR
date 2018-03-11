using Distributions

function bayesPR(genoTrain, phenoTrain, snpInfo, chrs, fixedRegSize, varGenotypic, varResidual, chainLength, burnIn, outputFreq, onScreen)
    SNPgroups = prepRegionData(snpInfo, chrs, genoTrain, fixedRegSize)
    fileControl(fixedRegSize)
    these2Keep = collect((burnIn+1):outputFreq:chainLength) #print these iterations
    nRegions    = length(SNPgroups)
    println("number of regions: ", nRegions)
    dfEffectVar = 4.0
    dfRes       = 4.0
    X           = convert(Array{Float64}, genoTrain)
    println("X is this size", size(X))
    y           = convert(Array{Float64}, phenoTrain)
    println("y is this size", size(y))
    nRecords , nMarkers   = size(y,1), size(X,2)
    p           = mean(X,1)/2.0
    sum2pq      = sum(2*(1-p).*p)
    varBeta         = fill(varGenotypic/sum2pq, nRegions)
    scaleVar        = varBeta[1]*(dfEffectVar-2.0)/dfEffectVar
    νS_β            = scaleVar*dfEffectVar
    df_β            = dfEffectVar
    scaleRes        = varResidual*(dfRes-2.0)/dfRes
    νS_e            = scaleRes*dfRes
    df_e            = dfRes
    tempBetaVec     = zeros(Float64,nMarkers) #initial values as "0"
    μ               = mean(y)
    X              .-= ones(Float64,nRecords)*2p
    xpx             = diag(X'X)
    ycorr           = y .- μ
    #MCMC starts here
    for iter in 1:chainLength
        #sample residual variance
        varE = sampleVarE(νS_e,ycorr,df_e,nRecords)
        #sample intercept
        ycorr    .+= μ
        rhs      = sum(ycorr)
        invLhs   = 1.0/nRecords
        meanMu   = rhs*invLhs
        μ        = rand(Normal(meanMu,sqrt(invLhs*varE)))
        ycorr    .-= μ
        for r in 1:nRegions
            theseLoci = SNPgroups[r]
            regionSize = length(theseLoci)
            λ_r = varE/varBeta[r]
            for l in theseLoci
                BLAS.axpy!(tempBetaVec[l], X[:,l], ycorr)
                rhs = X[:,l]'*ycorr
                lhs = xpx[l] + λ_r
                meanBeta = lhs\rhs                                     #I can use invLhs = 1.0/lhs as lhs is a scalar
                tempBetaVec[l] = sampleBeta(meanBeta, lhs, varE)       #then I can use invLhs*varE
                BLAS.axpy!(-1*tempBetaVec[l], X[:,l], ycorr)
            end
            varBeta[r] = sampleVarBeta(νS_β,tempBetaVec[theseLoci],df_β,regionSize)
        end
        outputControl(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    end
#    betaFromFile =  readcsv(pwd()"/betaOut",header=false)
#    print("read $(size(betaFromFile,1)) samples for $(size(betaFromFile,2)) markers from betaOut \n")
#    bayesRegOut = mean(betaFromFile,1)
#    @printf("acc %.6f \n", cor(y,X*bayesRegOut')[1])
end

function prepRegionData(snpInfo,chrs,genoTrain,fixedRegSize)
    println("prepRegionData started")
    accRegion = 0
    accRegionVec = [0]
    SNPgroups = []
#    mapData = readtable(pwd()"/$mapFile", header=false)
    ###only for our map file
    mapData = readtable("$snpInfo", header=false, separator=' ')
    headMap = [:row, :snpID, :snpOrder ,:chrID, :pos]
    rename!(mapData , names(mapData), headMap)
    println("line73")
    print(head(mapData))
    println("line73")
    ###
    mapData[:chrID] .<= chrs
    totLoci = size(genoTrain,2)
    snpInfoFinal = DataFrame(Any, 0, 3)
    for c in 1:chrs
        thisChr = mapData[mapData[:chrID] .== c,:]
        totLociChr = size(thisChr,1)
        TotRegions = ceil(Int,totLociChr/fixedRegSize)
        accRegion += TotRegions
        push!(accRegionVec, accRegion)
        tempGroups = sort(repmat(collect(accRegionVec[c]+1:accRegionVec[c+1]),fixedRegSize))
#        tempGroups = sort(repmat(collect(1:TotRegions),fixedRegSize))
        snpInfo = DataFrame(Any, length(tempGroups), 3)
        snpInfo[1:totLociChr,1] = collect(1:totLociChr)
        snpInfo[1:totLociChr,2] = thisChr[:snpID]
        snpInfo[:,3] = tempGroups
        dropmissing!(snpInfo)
        snpInfoFinal = vcat(snpInfoFinal,snpInfo)
#        rename!(snpInfo, names(snpInfo), [:snpOrder, :snpID, :regionID])
        @printf("chr %.0f has %.0f groups \n", c, TotRegions)
#        println(counts(snpInfo[:,3]))
            println("line97")
        println(by(snpInfo, :x3, nrow)[:,2])
            println("line99")
    end
    print(snpInfoFinal)
    for g in 1:accRegion
        push!(SNPgroups,searchsorted(snpInfoFinal[:,3], g))
    end
    println("prepRegionData ended")
    return SNPgroups
end

function outputControl(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    println("outputControl started")
    if iter in these2Keep
        out0 = open(pwd()*"/muOut$fixedRegSize", "a")
        writecsv(out0, μ)
        close(out0) 
        out1 = open(pwd()*"/betaOut$fixedRegSize", "a")
        writecsv(out1, tempBetaVec')
        close(out1)
        out2 = open(pwd()*"/varBetaOut$fixedRegSize", "a")
        writecsv(out2, varBeta')
        close(out2)
        out3 = open(pwd()*"/varEOut$fixedRegSize", "a")
        writecsv(out3, varE)
        close(out3)    
        if onScreen==true
            varU = var(X*tempBetaVec)
            @printf("iter %s varU %.2f varE %.2f\n", iter, varU, varE)
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
    println("outputControl ended")
end

function fileControl(fixedRegSize)
    println("fileControl started")
    for f in ["muOut$fixedRegSize" "betaOut$fixedRegSize" "varBetaOut$fixedRegSize" "varEOut$fixedRegSize"]
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
    println("fileControl ended")
end

function sampleBeta(meanBeta, lhs, varE)
    return rand(Normal(meanBeta,sqrt(lhs\varE)))
end

function sampleVarBeta(νS_β,whichLoci,df_β,regionSize)
    return((νS_β + dot(whichLoci,whichLoci))/rand(Chisq(df_β + regionSize)))
end
function sampleVarE(νS_e,yCorVec,df_e,nRecords)
    return((νS_e + dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords)))
end
