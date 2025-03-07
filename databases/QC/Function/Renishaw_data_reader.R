
###Renishaw单光谱

wavelength_summary <- function(hyperspec,Name_group)
{
  meta <- as.data.frame(tstrsplit(gsub(".txt","",hyperspec$filename),"_"),col.names = Name_group)
  meta <- mutate(meta,Group = apply(dplyr::select(meta,starts_with("group")),1,function(i) paste(i,collapse = "_")))
  hyperspec@data <- cbind(meta,hyperspec@data)
  wavelength_summary <- data.frame(Group = hyperspec$Group[1],
                                   first_wavnumber = hyperspec@wavelength[1],
                                   length = length(hyperspec@wavelength))
  return(wavelength_summary)
}

renishaw <- function(folderpath,Name_group,save = T,outputpath)
{
  filenames <- list.files(folderpath,pattern = "*.txt")
  ###将所有Txt的第一行删除
  for (filename in filenames)
  {temp <- fread(filename,header = T,sep = "\t")
  fwrite(temp,file = filename,sep = "\t") }
  ###预读取所有数据并根据Wavenumber的统计结果分组
  wavelengths <- read.txt.Renishaw(paste(folderpath,filenames[1],sep = "/"),data = "xyspc")
  wavelength_summarize <- wavelength_summary(wavelengths,Name_group)
  for (filename in filenames[-1])
  {
    temp <- read.txt.Renishaw(paste(folderpath,filename,sep = "/"),data = "xyspc")
    temp_summarize <- wavelength_summary(temp,Name_group)
    wavelength_summarize <- rbind(wavelength_summarize,temp_summarize)
  }
  wavelength_summarize$group <- paste(wavelength_summarize$first_wavnumber,wavelength_summarize$length)
  if (save == T)
  {
    ###根据分组结果读取数据
    for (i in levels(factor(wavelength_summarize$group)))
    {
      group_i <- which(wavelength_summarize$group == i)
      filenames_i <- filenames[group_i]
      final_dataframe <- read.txt.Renishaw(paste(folderpath,filenames_i[1],sep = "/"),data = "xyspc")
      for (filename in filenames_i)
      {
        temp <- read.txt.Renishaw(paste(folderpath,filename,sep = "/"),data = "xyspc")
        final_dataframe <- rbind(final_dataframe,temp)
      }
      data_hyperSpec <- final_dataframe
      data_hyperSpec$filename <- basename(data_hyperSpec$filename)
      data_hyperSpec$Number <- c(1:nrow(data_hyperSpec))
      meta <- as.data.frame(tstrsplit(gsub(".txt","",basename(data_hyperSpec$filename)),"_"),col.names = Name_group)
      meta <- mutate(meta,Group = apply(dplyr::select(meta,starts_with("group")),1,function(i) paste(i,collapse = "_")))
      data_hyperSpec@data <- cbind(meta,data_hyperSpec@data)
      single_cell(outputpath,"all_single_cell",data_hyperSpec)
    }
  }
  return(wavelength_summarize)
}

