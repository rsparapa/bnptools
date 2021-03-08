hbart_save = function(
                      ifilename, ##filename
                         ix, ##training points
                         im,  ##number of trees in mean model
                         imh, ##number of trees in variance model
                         ind, ##number of draws saved from the posterior
                         ixicuts, ##variable cutpoints for each predictor
                         ifit ##saved fitted model object returned 
                         ) {
    res=.Call("cpsambrt_save",
              ifilename, ##filename
              ix, ##training points
              im,  ##number of trees in mean model
              imh, ##number of trees in variance model
              ind, ##number of draws saved from the posterior
              ixicuts, ##variable cutpoints for each predictor
              ifit, ##saved fitted model object returned 
              PACKAGE='hbart')
    return(res)
}
