#!/usr/bin/env python
"""
Thomas:

"""

########################################
# Class
########################################

class TreeNode:

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.children = []

    def Print(self):

        print 'Class TreeNode: {0}'.format( self.split_name )
        print '    n = {0}'.format( self.n )

        if len(self.children) > 0:
            child_names = [ child.name for child in self.children ]
            print '    children = {0}'.format( child_names )
        else:
            print '    children = []'


########################################
# Functions
########################################

def Get_child( *input_args ):

    if len( input_args ) == 0:
        print "Format: Get_number( <TreeNode>, <1/0>, <1/0>, ... )"
        return 0

    t = input_args[0]

    for i in input_args[1:]:

        if len( t.children ) == 0:
            print 'Too many input arguments - Tree does not have that many children'
            return 0

        t = t.children[i]

    return t


def Print_Tree( t ):

    print '\nPrinting match tree {0}'.format( t.name )
    print '===================='
    print '  Total count         = {0}'.format( t.n )
    print '  Total passed to MEM = {0}\n'.format( Get_child(t,1,1,1).n )

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        'Jet to Qrk',
        '',
        'Subj to Qrk',
        '',
        'Jet to Subj',
        '',
        )

    print '|-------------+----------+-------------+----------+-------------+----------|'

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        'Success',
        str( Get_child(t,1).n ),
        'Success',
        str( Get_child(t,1,1).n ),
        'Success',
        str( Get_child(t,1,1,1).n ),
        )


    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        '',
        '',
        '',
        '',
        'Failed',
        str( Get_child(t,1,1,0).n ),
        )

    print '|             |          |-------------+----------+-------------+----------|'

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        '',
        '',
        'Failed',
        str( Get_child(t,1,0).n ),
        'Success',
        str( Get_child(t,1,0,1).n ),
        )

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        '',
        '',
        '',
        '',
        'Failed',
        str( Get_child(t,1,0,0).n ),
        )

    print '|-------------+----------+-------------+----------+-------------+----------|'

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        'Failed',
        str( Get_child(t,0).n ),
        'Success',
        str( Get_child(t,0,1).n ),
        'Success',
        str( Get_child(t,0,1,1).n ),
        )


    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        '',
        '',
        '',
        '',
        'Failed',
        str( Get_child(t,0,1,0).n ),
        )

    print '|             |          |-------------+----------+-------------+----------|'

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        '',
        '',
        'Failed',
        str( Get_child(t,0,0).n ),
        'Success',
        str( Get_child(t,0,0,1).n ),
        )

    print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
        '',
        '',
        '',
        '',
        'Failed',
        str( Get_child(t,0,0,0).n ),
        )

    print '----------------------------------------------------------------------------'


########################################
# Main
########################################

def main():

    t = TreeNode( 'b' )

    # 1st level
    t.children = [ TreeNode('split1_fail'), TreeNode('split1_success') ]

    # 2nd level
    b = Get_child( t, 0 )
    b.children = [ TreeNode('split2_fail'), TreeNode('split2_success') ]

    b = Get_child( t, 1 )
    b.children = [ TreeNode('split2_fail'), TreeNode('split2_success') ]

    # 3rd level
    b = Get_child( t, 0, 0 )
    b.children = [ TreeNode('split3_fail'), TreeNode('split3_success') ]

    b = Get_child( t, 0, 1 )
    b.children = [ TreeNode('split3_fail'), TreeNode('split3_success') ]

    b = Get_child( t, 1, 0 )
    b.children = [ TreeNode('split3_fail'), TreeNode('split3_success') ]

    b = Get_child( t, 1, 1 )
    b.children = [ TreeNode('split3_fail'), TreeNode('split3_success') ]

    Print_Tree( t )




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
