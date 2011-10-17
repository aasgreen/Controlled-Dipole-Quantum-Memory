#This script will generate a eps file from the psfrag generated matlab script. 
#I will extend this to work with variable names, but for now, just make sure that the matlabfrag generated files are called test1

latex plot_cmd.tex

#this will create the .dvi file. Now, we need to convert this to a ps

dvips -o fig.ps plot_cmd.dvi

ps2epsi fig.ps fig.epsi

#Strip preview part from fig.epsi

perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/' <fig.epsi > fig.eps
