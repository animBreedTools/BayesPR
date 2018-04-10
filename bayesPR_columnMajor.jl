using DataFrames
using Distributions

function bayesPR(genoTrain, phenoTrain, snpInfo, chrs, fixedRegSize, varGenotypic, varResidual, chainLength, burnIn, outputFreq, onScreen)
    SNPgroups, genoX = prepRegionData(snpInfo, chrs, genoTrain, fixedRegSize)
    these2Keep = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations
    nRegions    = length(SNPgroups)
    println("number of regions: ", nRegions)
    dfEffectVar = 4.0
    dfRes       = 4.0
    X           = convert(Array{Float64}, genoX[:,2:end])
    println("X is this size", size(X))
    y           = convert(Array{Float64}, phenoTrain)
    println("y is this size", size(y))
    nTraits, nRecords , nMarkers   = size(y,2), size(y,1), size(X,2)
    fileControlSt(fixedRegSize)
    p           = mean(X,1)./2.0
    sum2pq      = sum(2*(1-p).*p)
    if varGenotypic==0.0
        varBeta      = fill(0.0005, nRegions)
        else varBeta = fill(varGenotypic/sum2pq, nRegions)
    end
    if varResidual==0.0
        varResidual  = 0.0005
    end
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
        outputControlSt(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    end
#    betaFromFile =  readcsv(pwd()"/betaOut",header=false)
#    print("read $(size(betaFromFile,1)) samples for $(size(betaFromFile,2)) markers from betaOut \n")
#    bayesRegOut = mean(betaFromFile,1)
#    @printf("acc %.6f \n", cor(y,X*bayesRegOut')[1])
end

function mtBayesPR(genoTrain, phenoTrain, snpInfo, chrs, fixedRegSize, varGenotypic, varResidual, chainLength, burnIn, outputFreq, onScreen)
    SNPgroups, genoX = prepRegionData(snpInfo, chrs, genoTrain, fixedRegSize)
    these2Keep = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations
    nRegions    = length(SNPgroups)
    println("number of regions: ", nRegions)
    dfEffect    = 4.0
    dfRes       = 4.0
    X           = convert(Array{Float64}, genoX[:,2:end])  #first colum is ID
    println("X is this size", size(X))
    Y           = convert(Array{Float64}, phenoTrain)
    println("Y is this size", size(Y))
    nTraits, nRecords , nMarkers   = size(Y,2), size(Y,1), size(X,2)
    fileControl(nTraits,fixedRegSize)
    p           = mean(X,1)./2.0
    sum2pq      = sum(2*(1-p).*p)
    if varGenotypic==0.0
        covBeta       = fill([0.003 0;0 0.003],nRegions)
        else covBeta  = fill(varGenotypic/sum2pq,nRegions)
    end
    if varResidual==0.0
        varResidual   = [0.003 0;0 0.003]
    end
     #priors#
    dfβ         = dfEffect + nTraits
    dfR         = dfRes + nTraits
    Vb          = covBeta[1].*(dfEffect-nTraits-1)
    VR          = varResidual.*(dfR - nTraits - 1)
    #initial Beta values as "0"
    tempBetaMat     = zeros(Float64,nTraits,nMarkers)
    μ = mean(Y,1)    
    X              .-= ones(Float64,nRecords)*2p
    xpx             = diag(X'X)
    #allocate memory for mmeRun2
#    x1x2    = zeros(Float64,length(vec(Y)),nTraits)
#    x1x2pRi = zeros(Float64,nTraits,length(vec(Y)))
#    x1x2T   = Array{Array{Float32,2}}(nMarkers)
#    for locus in 1:nMarkers
#        x1x2T[locus] = [X[:,locus] zeros(nRecords);zeros(nRecords) X[:,locus]]
#    end

    ycorr1 = zeros(Float64,nTraits,length(Y[:,1]))
    ycorr2 = zeros(Float64,nTraits,length(Y[:,2]))
    
    ycorr1 = (Y[:,1] .- μ[1])
    ycorr2 = (Y[:,2] .- μ[2])
        
    for iter in 1:chainLength
        #sample residual var
        Rmat = sampleCovarE(dfR, nRecords, VR, ycorr1, ycorr2)
#        Ri = kron(inv(Rmat),eye(nRecords)) #for mmeRun2
        Ri = inv(Rmat)

        # sample intercept
        ycorr1 += μ[1]
        rhs1 = sum(ycorr1)
        invLhs1 = 1.0/nRecords
        mean1 = rhs1*invLhs1

        #    sample intercept
        ycorr2 += μ[2]
        rhs2 = sum(ycorr2)
        invLhs2 = 1.0/nRecords
        mean2 = rhs2*invLhs2
    
        μ[1] = rand(Normal(mean1,sqrt(invLhs1*Rmat[1,1])))
        μ[2] = rand(Normal(mean2,sqrt(invLhs2*Rmat[2,2])))
        
        ycorr1 -= μ[1]
        ycorr2 -= μ[2]
        
        for r in 1:nRegions
            theseLoci = SNPgroups[r]
            regionSize = length(theseLoci)
            invB = inv(covBeta[r])
            for locus in theseLoci
                BLAS.axpy!(tempBetaMat[1,locus], X[:,locus], ycorr1)
                BLAS.axpy!(tempBetaMat[2,locus], X[:,locus], ycorr2) 
                tempBetaMat[:,locus] = mmeRunFast(X[:,locus],Ri,locus,xpx,ycorr1,ycorr2,invB)
                BLAS.axpy!(-1*tempBetaMat[1,locus], X[:,locus], ycorr1)
                BLAS.axpy!(-1*tempBetaMat[2,locus], X[:,locus], ycorr2)
            end
            covBeta[r] = sampleCovBeta(dfβ,regionSize,Vb,tempBetaMat, theseLoci)
        end
        outputControl(nTraits,onScreen,iter,these2Keep,X,tempBetaMat,μ,covBeta,Rmat,fixedRegSize,nRegions)
    end
end

function prepRegionData(snpInfo,chrs,genoTrain,fixedRegSize)
    accRegion = 0
    accRegionVec = [0]
    SNPgroups = []
#    mapData = readtable(pwd()"/$mapFile", header=false)
    ###only for our map file
    mapData = readtable("$snpInfo", header=false, separator=' ')
    headMap = [:row, :snpID, :snpOrder ,:chrID, :pos]
    rename!(mapData , names(mapData), headMap)
    print(mapData[1:5,:])
    mapData[:snpID] = ["M$i" for i in 1:size(mapData,1)] #to convert original IDs like "HAPMAP43437-BTA-101873"
    print(mapData[1:10,:])
    ###
    mapData = mapData[mapData[:chrID] .<= chrs,:]
    # if first col in genoTrain is ID
    # I find cols that are in mapData (<chrs), and select those
    usedLoci = intersect(names(genoTrain),Symbol.(mapData[:snpID]))
    println("usedLoci",usedLoci)
    mapData = mapData[[find(usedLoci[i].==Symbol.(mapData[:snpID]))[] for i in 1:length(usedLoci)],:] #trim map data
    println("mapData",mapData)
    genoX = genoTrain[:,vcat(Symbol("ID"),usedLoci)]    #trim genoData
#     genoX = genoTrain[:,[1; [find(i -> i == j, names(genoTrain))[] for j in [Symbol(mapData[:snpID][i]) for i in 1:size(mapData,1)]]]]
    #genoX = genoTrain[:,[find(i -> i == j, names(genoTrain))[] for j in [Symbol(mapData[:snpID][i]) for i in 1:size(mapData,1)]]]
    totLoci = size(genoX[:,2:end],2) # first col is ID
    snpInfoFinal = DataFrame(Any, 0, 3)
    if fixedRegSize==99
        println("fixedRedSize $fixedRegSize")
        snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
        accRegion    = length(unique(mapData[:chrID]))
        elseif fixedRegSize==9999
            snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
            snpInfoFinal[:,:chrID]  .= 1
            accRegion    = 1
        else
        for c in 1:chrs
            thisChr = mapData[mapData[:chrID] .== c,:]
            totLociChr = size(thisChr,1)
            TotRegions = ceil(Int,totLociChr/fixedRegSize)
            accRegion += TotRegions
            push!(accRegionVec, accRegion)
            tempGroups = sort(repmat(collect(accRegionVec[c]+1:accRegionVec[c+1]),fixedRegSize))
#           tempGroups = sort(repmat(collect(1:TotRegions),fixedRegSize))
            snpInfo = DataFrame(Any, length(tempGroups), 3)
            snpInfo[1:totLociChr,1] = collect(1:totLociChr)
            snpInfo[1:totLociChr,2] = thisChr[:snpID]
            snpInfo[:,3] = tempGroups
            dropmissing!(snpInfo)
            snpInfoFinal = vcat(snpInfoFinal,snpInfo)
#           rename!(snpInfo, names(snpInfo), [:snpOrder, :snpID, :regionID])
            @printf("chr %.0f has %.0f groups \n", c, TotRegions)
#           println(counts(snpInfo[:,3]))
            println(by(snpInfo, :x3, nrow)[:,2])
        end
        end  #ends if control flow
    print(snpInfoFinal)
    writecsv("snpInfo",convert(Array,snpInfoFinal))
    for g in 1:accRegion
        push!(SNPgroups,searchsorted(snpInfoFinal[:,3], g))
    end
    return SNPgroups, genoX
end

function outputControlSt(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
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
end

function fileControlSt(fixedRegSize)
    for f in ["muOut$fixedRegSize" "betaOut$fixedRegSize" "varBetaOut$fixedRegSize" "varEOut$fixedRegSize"]
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
end

function outputControl(nTraits,onScreen,iter,these2Keep,X,tempBetaMat,μ,covBeta,Rmat,fixedRegSize,nRegions)
    if iter in these2Keep
        out0 = open(pwd()*"/muOut$fixedRegSize", "a")
        writecsv(out0, μ)
        close(out0)
        for t in 1:nTraits
            out1 = open(pwd()*"/beta"*"$t"*"Out$fixedRegSize", "a")
            writecsv(out1, tempBetaMat[t,:]')
            close(out1)
            out2 = open(pwd()*"/varBeta"*"$t"*"Out$fixedRegSize", "a")
            printThis = [vcat(covBeta[r]...)[t^2] for r in 1:nRegions]'
            writecsv(out2, printThis) #works only for bivariate
            close(out2)
        end
        outCov = open(pwd()*"/covBetaOut$fixedRegSize", "a")
        printThis = [vcat(covBeta[r]...)[3] for r in 1:nRegions]'
        writecsv(outCov, printThis) #works only for bivariate
        close(outCov)
        out3 = open(pwd()*"/varEOut$fixedRegSize", "a")
        writecsv(out3, vec(Rmat)')
        close(out3)    
        if onScreen==true
            varU = cov(X*tempBetaMat')
            println("iter $iter \nvarU: $varU \nvarE: $Rmat \n")
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
end

function fileControl(nTraits,fixedRegSize)
    files2Remove = ["muOut$fixedRegSize","varEOut$fixedRegSize","covBetaOut$fixedRegSize"]
    for t in 1:nTraits
        push!(files2Remove,"beta"*"$t"*"Out$fixedRegSize")
        push!(files2Remove,"varBeta"*"$t"*"Out$fixedRegSize")
    end
    for f in files2Remove
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
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
function sampleCovBeta(dfβ, regionSize, Vb , tempBetaMat, theseLoci)
    Sb = tempBetaMat[:,theseLoci]*tempBetaMat[:,theseLoci]'
    return rand(InverseWishart(dfβ + regionSize, Vb + Sb))
end
function sampleCovarE(dfR, nRecords, VR, ycorr1, ycorr2)
     Sr = [ycorr1 ycorr2]'* [ycorr1 ycorr2]
    return rand(InverseWishart(dfR + nRecords,VR + Sr))
end
function mmeRun1(X,locus,nRecords,ycorr1,ycorr2,Ri,invB)
    x1x2   = ([X[:,locus] zeros(nRecords);zeros(nRecords) X[:,locus]])
    x1x2pRi= x1x2'*Ri
    rhs    = x1x2pRi*[ycorr1;ycorr2]
    invLhs = inv(x1x2pRi*x1x2 + invB)
    meanBeta = invLhs*rhs    
    return rand(MvNormal(meanBeta,convert(Array,Symmetric(invLhs))))
end
#function mmeRun2(locus,x1x2,x1x2T,x1x2pRi,nRecords,ycorr1,ycorr2,Ri,invB) #memory optimized way
#    x1x2 .= x1x2T[locus]
#    x1x2pRi  .= x1x2'*Ri
#    rhs    = x1x2pRi*[ycorr1;ycorr2]
#    invLhs = inv(x1x2pRi*x1x2 + invB)
#    meanBeta = invLhs*rhs    
#    return rand(MvNormal(meanBeta,convert(Array,Symmetric(invLhs))))
#end
function mmeRunFast(x,Ri,locus,xpx,ycorr1,ycorr2,invB)
    r1 = x'*Ri[1]
    r2 = x'*Ri[2]
    r3 = x'*Ri[3]
    r4 = x'*Ri[4]

    rhs    = [r1 r2;r3 r4]*[ycorr1;ycorr2]
    invLhs = inv([xpx[locus] xpx[locus];
                  xpx[locus] xpx[locus]].*Ri + invB)
    
    meanBeta = invLhs*rhs    
    return rand(MvNormal(meanBeta,convert(Array,Symmetric(invLhs))))
end

