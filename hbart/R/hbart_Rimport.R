hbart_Rimport = function(
                         ix, ##training points
                         im,  ##number of trees in mean model
                         imh, ##number of trees in variance model
                         ind, ##number of draws saved from the posterior
                         ixicuts, ##variable cutpoints for each predictor
                         ipoststr ##saved posterior tree draws
                         )
{
    res=.Call("cpsambrt_Rimport",
        ix, ##training points
        im,  ##number of trees in mean model
        imh, ##number of trees in variance model
        ind, ##number of draws saved from the posterior
        ixicuts, ##variable cutpoints for each predictor
        ipoststr, ##saved posterior tree draws
        PACKAGE='hbart')
    return(res)
}
