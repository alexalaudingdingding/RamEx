###RamEx installer
###Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
###Updated at Sep. 22, 2020 
###Updated by Gongchao Jing
#!/bin/bash
##Users can change the default environment variables configuration file here
if [[ $SHELL = '/bin/zsh' ]];
then
        PATH_File=~/.zshrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.zsh_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi
else
        PATH_File=~/.bashrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.bash_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi

fi
RamEx_PATH=`pwd`
Sys_ver=`uname`
###Checking that environment variable of RamEx exists###
Check_old_ramex=`grep "export RamEx"  $PATH_File|awk -F '=' '{print $1}'`
Check_old_path=`grep "RamEx/bin"  $PATH_File |sed 's/\(.\).*/\1/' |awk '{if($1!="#"){print "Ture";}}'`
Add_Part="####DisabledbyRamEx####"
echo "**RamEx Installation**"
echo "**version 1.0**"

###Build source code for src package###
if [ -f "Makefile" ]
   then
       echo -e "\n**RamEx src package**"
       make
       echo -e "\n**Build Complete**"
else
   echo -e "\n**RamEx bin package**"
fi
###Configure environment variables###

if [ "$Check_old_ramex" != "" ]
   then
      Checking=`grep ^export\ RamEX  $PATH_File|awk -F '=' '{print $2}'`
      if [ "$Checking" != "$RamEx" ]
         then
         if [ "$Sys_ver" = "Darwin" ]
            then
            sed -i "" "s/^export\ RamEx/$Add_Part\ &/g" $PATH_File
            sed -i "" -e "`grep -n "$Add_Part" $PATH_File | cut -d ":" -f 1 | head -1` a\ 
export\ RamEX=$RamEx
" $PATH_File
         else
             sed -i "s/^export\ RamEx/$Add_Part\ &/g" $PATH_File
             sed -i "/$Add_Part\ export\ RamEx/a export\ RamEx=$RamEx" $PATH_File
         fi
     fi    
elif [ "$Check_old_ramex" = "" ]
    then
      echo "export RamEx="${RamEx} >> $PATH_File
fi
if [ "$Check_old_path" = "" ]
    then
      echo "export PATH=\$PATH:\$RamEx/bin" >> $PATH_File
fi
###Source the environment variable file###
source $PATH_File
echo -e "\n**Environment Variables Configuration Complete**"
###Configurate the R packages###
echo -e ""
Rscript $RamEx_PATH/Rscript/config.R
###End
echo -e "\n**RamEx Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"



