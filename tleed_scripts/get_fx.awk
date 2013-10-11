awk '
#/TRUTH/{print $0}
/fx =/{lin7=NR}
        {if (NR ==lin7+2 && lin7 > 0)
        print $1}
' tmp.log 
