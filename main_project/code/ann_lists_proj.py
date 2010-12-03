import sys
import os;

for arg in sys.argv[1:len(sys.argv)]:
	for x in range(0,5):
		xstr = str(x)
		aa = 'SMM/' + arg + '/c00' + xstr + ' | grep -v "#" >'
		os.system( 'src/seq2inp '+aa+' sp_c00' + xstr )
		os.system( 'src/seq2inp -bl '+aa+' b_c00' + xstr )

		bb = 'SMM/' + arg + '/c00' + xstr + ' | grep -v "#" >'
		os.system( 'src/seq2inp '+bb+' sp_f00' + xstr )
		os.system( 'src/seq2inp -bl '+bb+' b_f00' + xstr )

#		for y in range(x+1,5):
#			ystr = str(y)
#			for m in range(0,5):
#				if m != x and m != y:
#					os.system( 'src/seq2inp '+aa+'> sp_f00' + xstr + '_' + ystr )
#					os.system( 'src/seq2inp -bl '+aa+'> b_f00' + xstr + '_' + ystr )

	for x in range(0,5):
		xstr = str(x)
		for nh in [1, 2, 10]:
			nhstr = str(nh)
			os.system( 'src/nnbackprop -nh '+nhstr+' -tf sp_c00'+xstr+' -syn '+nhstr+'_sp_syn_dat_'+xstr+' sp_f00'+xstr+' > /dev/null' )
			os.system( 'src/nnbackprop -nh '+nhstr+' -tf b_c00'+xstr+' -syn '+nhstr+'_b_syn_dat_'+xstr+' b_f00'+xstr+' > /dev/null' )
			
	for x in range(0,5):
		xstr = str(x)
		for nh in [1, 2, 10]:
			nhstr = str(nh)
			for u in range(0,5):
				if u != x:
					ustr = str(u)	
					os.system( 'echo '+nhstr+'_sp_syn_dat_'+ustr+' >> list_'+nhstr+'_sp'+xstr )
					os.system( 'echo '+nhstr+'_b_syn_dat_'+ustr+' >> list_'+nhstr+'_b'+xstr )

	for x in range(0,5):
		xstr = str(x)
		for nh in [1, 2, 10]:
			nhstr = str(nh)
			print( 'NN '+nhstr+'HL\t'+arg+'\t'+xstr )
			run = 'src/nnforward list_'+nhstr+'_sp'+xstr+' sp_c00'+xstr+' | grep -v "#" | gawk \'{print $1,$3}\' | src/xycorr | gawk \'{OFS = "\t" ; print $5,$7}\''
			os.system( run )
			print( 'NN '+nhstr+'HL Blosum\t'+arg+'\t'+xstr )
			run = 'src/nnforward list_'+nhstr+'_b'+xstr+' b_c00'+xstr+' | grep -v "#" | gawk \'{print $1,$3}\' | src/xycorr | gawk \'{OFS = "\t" ; print $5,$7}\''
			os.system( run )
