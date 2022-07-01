#' @title br_harmony
#'
#' @description use Harmony to align the test data and training data onto a PCA
#'
#' @param train.xx A data matrix where each row is a cell and each column is a gene in training data
#' @param test A data matrix where each row is a cell and each column is a gene  in test data
#'
#' @return A list containing aligned PCs matrix for training data and test data.
#'
#' @importFrom harmony HarmonyMatrix
#' @importFrom Seurat CreateSeuratObject RunPCA ScaleData
#' @noRd
#'
br_harmony=function(train.xx,test){

  meta_data=data.frame(dataset=c(rep("train",nrow(train.xx)),rep("test",nrow(test))))
  #pca
  obj.pc=CreateSeuratObject(counts = cbind(t(train.xx),t(test)))
  obj.pc=FindVariableFeatures(obj.pc)
  obj.pc=ScaleData(obj.pc)
  obj.pc=RunPCA(obj.pc,npcs = 20)

  pc.matrix= obj.pc@reductions$pca@cell.embeddings

  my_harmony_embeddings = HarmonyMatrix(
    data_mat=pc.matrix,
    meta_data=meta_data,
    vars_use="dataset",
    do_pca=FALSE
  )

  train.xx_new=my_harmony_embeddings[meta_data$dataset=="train",]
  test_new=my_harmony_embeddings[meta_data$dataset=="test",]

  oup=list(train=train.xx_new,
           test=test_new)
  return(oup)
}
