import sys
from optparse import OptionParser, IndentedHelpFormatter

def run(n):
    '''takes in some number n, returns nth number of the fibonacci sequence. Uses
        a recursive function to accomplish this'''
    n=int(n)
    if n<2:
        return n
    return run(n-1)+run(n-2)




def main(args):
    
	parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
	parser.set_description(main.__doc__)

    	parser.add_option('--iterations', '-n',
    action="store", default=None,
    help="n iterations of Fibonacci Sequence", )

	#parse option
	(options, args) = parser.parse_args(args=args[1:])

	global Options
	Options = options

        x=run(Options.iterations);
        print x

if __name__ == "__main__" : main(sys.argv)
