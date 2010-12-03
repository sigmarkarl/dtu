import sys
import os;

for arg in sys.argv[1:len(sys.argv)]:
	for x in range(0,5):
		xstr = str(x)
		aa = 'SMM/' + arg + '/c00' + xstr + ' | grep -v "#" >'
		os.system( 'src/seq2inp '+aa+' sp_c00' + xstr )
		os.system( 'src/seq2inp -bl '+aa+' b_c00' + xstr )

#		bb = 'SMM/' + arg + '/c00' + xstr + ' | grep -v "#" >'
#		os.system( 'src/seq2inp '+bb+' sp_f00' + xstr )
#		os.system( 'src/seq2inp -bl '+bb+' b_f00' + xstr )

		for y in range(x+1,5):
			ystr = str(y)
			for m in range(0,5):
				if m != x and m != y:
					mstr = str(m)
					mm = 'SMM/' + arg + '/c00' + mstr + ' | grep -v "#" >'
					os.system( 'src/seq2inp '+mm+'> sp_f00' + xstr + '_' + ystr )
					os.system( 'src/seq2inp -bl '+mm+'> b_f00' + xstr + '_' + ystr )

	for x in range(0,5):
		xstr = str(x)
		for nh in [1, 2, 10]:
			nhstr = str(nh)
			for y in range(0,5):
				ystr = str(y)
				for m in range(0,5):
					if m != x and m != y and x != y:
						xystr = str(min(x,y)) + '_' + str(max(x,y))
						nxystr = xstr + '_' + ystr
						os.system( 'src/nnbackprop -nh '+nhstr+' -tf sp_c00'+ystr+' -syn '+nhstr+'_sp_syn_dat_'+nxystr+' sp_f00'+xystr+' > /dev/null' )
						os.system( 'src/nnbackprop -nh '+nhstr+' -tf b_c00'+ystr+' -syn '+nhstr+'_b_syn_dat_'+nxystr+' b_f00'+xystr+' > /dev/null' )
			
	for x in range(0,5):
		xstr = str(x)
		for nh in [1, 2, 10]:
			nhstr = str(nh)
			for u in range(0,5):
				if u != x:
					ustr = str(u)
					os.system( 'echo '+nhstr+'_sp_syn_dat_'+xstr+'_'+ustr+' >> list_'+nhstr+'_sp'+xstr )
					os.system( 'echo '+nhstr+'_b_syn_dat_'+xstr+'_'+ustr+' >> list_'+nhstr+'_b'+xstr )

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

#os.system( "rm list*" )
#os.system( "rm sp_*" )
#os.system( "rm b_*" )
