#
# Script to make the wrfncxnj tarball
#
revision=$(svn info | grep Revision | awk '{print $2}')
mkdir wrfncxnj-0.1_r${revision}
cp README wrfncxnj-0.1_r${revision}
cp wrfncxnj.table_public wrfncxnj-0.1_r${revision}/wrfncxnj.table
for f in wrfncxnj.py wrfncxnj_base.py wrfncxnj_fun.py wrfncxnj_cli.py
  do sed -e "s/@progname@/$f/" gnu.lic | cat - $f > wrfncxnj-0.1_r${revision}/${f}
done

tar czvf wrfncxnj-0.1_r${revision}.tar.gz wrfncxnj-0.1_r${revision}/* 
