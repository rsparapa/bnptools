hbart_Rexport = function(
                         ix, ##training points
                         im,  ##number of trees in mean model
                         imh, ##number of trees in variance model
                         ind, ##number of draws saved from the posterior
                         ixicuts, ##variable cutpoints for each predictor
                         ifit ##saved fitted model object returned from cpsambrt
                         )
{
    res=.Call("cpsambrt_Rexport",
              ix, ##training points
              im,  ##number of trees in mean model
              imh, ##number of trees in variance model
              ind, ##number of draws saved from the posterior
              ixicuts, ##variable cutpoints for each predictor
              ifit, ##saved fitted model object returned from cpsambrt
              PACKAGE='hbart')
    return(res)
}
