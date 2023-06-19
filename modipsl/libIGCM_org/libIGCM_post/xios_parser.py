#! /usr/bin/env python
# coding: utf-8

#**************************************************************
# Author: Sebastien Denvil
# Contact: Sebastien.Denvil__at__ipsl.jussieu.fr
# $Revision:: 1552                                    $ Revision of last commit
# $Author:: aclsce                                    $ Author of last commit
# $Date:: 2022-04-08 09:28:51 +0000 (ج, 08 أبري�#$ Date of last commit
# IPSL (2006)
#  This software is governed by the CeCILL licence see libIGCM/libIGCM_CeCILL.LIC
#
#**************************************************************

import os, pwd, sys, traceback, argparse
import xml.etree.ElementTree as ET
#import readline, rlcompleter
#readline.parse_and_bind("tab: complete")

indent = 0
currentDepth = 0
ignoreElems = []
fromField=[]
fromFile=[]

def dump(args):
    """Dump XIOS xmls files."""
    # Read and dump xios_def_xml
    for inputFile in args.file:        
        # Read the file_def_xml
        print '\nReading %s \n|' % (inputFile)
        try:
            tree = ET.parse(inputFile)
        except:
            print "Parse error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return
        root=tree.getroot()
        # Call the recursive print
        printRecur(root)

def tsquery(args):
    """query timeseries related parameters from an XIOS xml file."""
    if args.verbosity >= 1: print 'Reading timeseries_def_xml=',args.file[0]
    try:
        tree = ET.parse(args.file[0])
    except:
        print "Parse error. Please fix so that it can be parsed."
        traceback.print_exc(file=sys.stdout)
        return
    root=tree.getroot()
    if args.verbosity >= 3: print root.tag, root.attrib
    findTimeSeries(root)
        
def printRecur(root):
    """Recursively prints the tree."""
    global indent
    global currentDepth
    if root.tag in ignoreElems:
        return
    print ' '*indent + '|--> %s: %s' % (root.tag, root.attrib)
    currentDepth += 1
    indent += 4
    if currentDepth <= args.depth or args.depth == None:
        for elem in list(root):
            printRecur(elem)
    currentDepth -= 1
    indent -= 4

def findTimeSeries(root):
    """Recursively find and list field tag with "timeseries", "id", "output_freq" and enable=.TRUE. attribute."""
    if root.tag == 'file' and root.attrib.get('timeseries'):
        if root.attrib.get('enabled') == '.TRUE.':
            print 'output_freq=%-5s,id=%s' % (root.attrib.get('output_freq'), root.attrib.get('id'))
    elif root.tag == 'file' and root.attrib.get('uuid_name'):
        print 'output_freq=%-5s,name=%85s' % (root.attrib.get('output_freq'), root.attrib.get('name'))
    else:
        for elem in list(root):
            findTimeSeries(elem)
    
def findField(root):
    """Recursively find and list field tag with "id" or "field_ref" attribute."""
    global fromField
    global fromFile
    if root.tag in ignoreElems:
        return    
    if root.tag == 'field' and root.attrib.get('id'):
        fromField.append(root.attrib.get('id'))
    elif root.tag == 'field' and root.attrib.get('field_ref'):
        fromFile.append(root.attrib.get('field_ref'))
    else:
        for elem in list(root):
            findField(elem)

def findFieldToRemove(root, fieldToRemove):
    """Recursively find tag having a field_ref in fieldToRemove."""
    if args.verbosity >= 3 and root.tag == 'file' and root.attrib.get('id'):
        print '\nFIELDS FROM FILE_DEF with id', root.attrib.get('id')
    for field in root.findall('field'):
        if args.verbosity >= 3: print 'field_ref=', field.attrib.get('field_ref')
        if field.attrib.get('field_ref') in fieldToRemove:
            if args.correction:
                if args.verbosity >= 2: print 'removing ', field.attrib.get('field_ref')
                root.remove(field)
            else:
                if args.verbosity >= 2: print 'To be removed ', field.attrib.get('field_ref')
    for elem in list(root):
        findFieldToRemove(elem, fieldToRemove)
        
def check(args):
    global fromFile
    exitCode=0
    # Read the field_def_xml
    if args.verbosity >= 1: print '\nReading field_def_xml=',args.field[0]
    try:
        tree = ET.parse(args.field[0])
    except:
        print 'Parse error with %s. Please fix so that it can be parsed.' % (args.field[0])
        traceback.print_exc(file=sys.stdout)
        return
    root=tree.getroot()
    if args.verbosity >= 3: print root.tag, root.attrib, '\n'
    # Build a list of field.id from field_def
    findField(root)

    # Loop over file_def files
    for inputFile in args.file:        
        # Read the file_def_xml
        if args.verbosity >= 1: print '\nReading file_def_xml=',inputFile
        try:
            tree = ET.parse(inputFile)
        except:
            print "Parse error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return
        root=tree.getroot()
        fromFile=[]
        if args.verbosity >= 3: print root.tag, root.attrib, '\n'

        # Build a list of field_ref from file_def
        findField(root)
        #print '4. fromFile=', fromFile
        
        # Compare the two lists. fromField must be a superset of fromFile.
        if set(fromField).issuperset(set(fromFile)):
            if args.verbosity >= 1: print '\nALL GOOD with %s' % (inputFile)
            if args.verbosity >= 3: print 'fromField=', fromField
            if args.verbosity >= 3: print 'fromFile=', fromFile
        else:
            if args.verbosity >= 1: print '\nTROUBLE AHEAD with %s' % (inputFile)
            if args.verbosity >= 3: print ', '.join(sorted(list(set(fromFile)-set(fromField))))
            # Identify fields in fromFile but not in fromField
            fieldToRemove=list(set(fromFile)-set(fromField))
            if args.verbosity >= 3: print 'fieldToRemove=', fieldToRemove
            # 
            # And now locate and remove them if the modify command has been called
            findFieldToRemove(root, fieldToRemove)
            # Final steps
            if args.correction: tree.write('modified.'+inputFile)
            if not len(fieldToRemove) == 0 and not args.correction:
                exitCode=1
    # The end
    sys.exit(exitCode)

def showtime(args):
    """
    prints table of formatted text format options
    """
    for style in xrange(6):
        for fg in xrange(30,36):
            s1 = ''
            for bg in xrange(40,46):
                format = ';'.join([str(style), str(fg), str(bg)])
                s1 += '\x1b[%sm %s \x1b[0m' % (format, pwd.getpwuid(os.getuid())[4]+' is on fire')
            print s1
        print '\n'
   
def modifydr2xmlIPSL(args):
    """Recursively find tree element with file tag and change the path prefix of the name attribute."""

    # Loop over file_def files
    for inputFile in args.file:
        # Read the file_def_xml
        if args.verbosity >= 1: print '\nReading dr2xml_def_xml=',inputFile
        try:
            tree = ET.parse(inputFile)
        except:
            print "Parse error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return
        root=tree.getroot()
        if args.verbosity >= 3: print root.tag, root.attrib, '\n'

        # Change file name prefix to point where we want
        for elem in tree.iter(tag='file'):
            if not elem.attrib.get('mode') == 'read':
                name=elem.attrib.get('name')
                if args.JobName is not None:
                    name = name.replace('JobName', args.JobName[0]).replace('ModelName', args.ModelName[0])
                    elem.set('name', os.path.normpath(args.newPath[0]) + '/' + os.path.normpath(name))
                else:
                    elem.set('name', os.path.normpath(args.newPath[0]) + '/' + os.path.normpath(name))
                if args.verbosity >= 2: print elem.tag, elem.attrib
                if args.splitLastDate is not None: elem.set('split_last_date', args.splitLastDate[0] + '-01-01 00:00:00')
                if args.ModelVersion is not None:
                    for childElem in elem.iter(tag='variable'):
                        if args.verbosity >= 2: print childElem.tag, childElem.attrib
                        if childElem.attrib.get('name') is not None:
                            if childElem.attrib.get('name') == 'model_version':
                               childElem.text = args.ModelVersion[0]
                if args.ModelName is not None:
                    for childElem in elem.iter(tag='variable'):
                        if args.verbosity >= 2: print childElem.tag, childElem.attrib
                        if childElem.attrib.get('name') is not None:
                            if childElem.attrib.get('name') == 'source_id':
                               childElem.text = args.ModelName[0]
                if args.ExperimentName is not None:
                    for childElem in elem.iter(tag='variable'):
                        if args.verbosity >= 2: print childElem.tag, childElem.attrib
                        if childElem.attrib.get('name') is not None:
                            if childElem.attrib.get('name') == 'experiment':
                               childElem.text = args.ExperimentName[0]


        # Write out the results
        if args.verbosity >= 1: print '\nWriting dr2xml_def_xml=','modified.' + inputFile
        try:
            tree.write ('modified.' + inputFile)
        except:
            print "Write error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return

def modifyPath(args):
    """Recursively find tree element with file tag and change the path prefix of the name attribute."""
    # Loop over file_def files
    for inputFile in args.file:
        # Read the file_def_xml
        if args.verbosity >= 1: print '\nReading dr2xml_def_xml=',inputFile
        try:
            tree = ET.parse(inputFile)
        except:
            print "Parse error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return
        root=tree.getroot()
        if args.verbosity >= 3: print root.tag, root.attrib, '\n'
        # Change file name prefix to point where we want
        for elem in tree.iter(tag='file'):
            if not elem.attrib.get('mode') == 'read':
                name=elem.attrib.get('name')
                if args.filePrefix is not None:
                    elem.set('name', os.path.normpath(args.newPath[0]) + '/' + args.filePrefix[0] + '_' + os.path.normpath(name))
                else:
                    elem.set('name', os.path.normpath(args.newPath[0]) + '/' + os.path.normpath(name))
                if args.verbosity >= 2: print elem.tag, elem.attrib
                   
        # Write out the results
        if args.verbosity >= 1: print '\nWriting dr2xml_def_xml=','modified.' + inputFile
        try:
            tree.write ('modified.' + inputFile)
        except:
            print "Write error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return

def modifyDR2XML(args):
    """Recursively find tree element with file tag and change the path prefix of the name attribute."""
    # Need to fine tune their freq_op and freq_offset
    cospTuple=('CMIP6_albisccp', 'CMIP6_cfadDbze94', 'CMIP6_cfadLidarsr532', 'CMIP6_clcalipso2', 'CMIP6_clcalipsoice', 'CMIP6_climodis', 'CMIP6_clisccp', 'CMIP6_clmisr', 'CMIP6_cltmodis', 'CMIP6_clwmodis', 'CMIP6_jpdftaureicemodis', 'CMIP6_jpdftaureliqmodis', 'CMIP6_pctisccp', 'CMIP6_parasolRefl', 'CMIP6_clcalipso', 'CMIP6_clhcalipso', 'CMIP6_cllcalipso', 'CMIP6_clmcalipso', 'CMIP6_cltcalipso', 'CMIP6_cltisccp')
    # Dont touch those when fine tuning cfsites field freq_op
    excluded_field=('CMIP6_ap', 'CMIP6_b', 'CMIP6_b_bnds', 'CMIP6_ap_bnds')
    
    # Loop over file_def files
    for inputFile in args.file:
        write=False
        # Read the file_def_xml
        if args.verbosity >= 1: print '\nReading dr2xml_def_xml=',inputFile
        try:
            tree = ET.parse(inputFile)
        except:
            print "Parse error. Please fix so that it can be parsed."
            traceback.print_exc(file=sys.stdout)
            return
        root=tree.getroot()
        if args.verbosity >= 3: print root.tag, root.attrib, '\n'

        # Change some elements to overrule dr2xml
        for elem in tree.iter(tag='field'):
            if elem.attrib.get('id') is not None:
                if elem.attrib.get('id') == 'CMIP6_siconc_Scaltypesi':
                    write=True
                    elem.set('freq_op', '2ts')
                    elem.set('freq_offset', '0ts')
                    if args.verbosity >= 2: print elem.tag, elem.attrib
                elif (elem.attrib.get('id').startswith(cospTuple)):
                    write=True
                    elem.set('freq_op', '3h')
                    elem.set('freq_offset', '0ts')
                    if args.verbosity >= 2: print elem.tag, elem.attrib
        # Change elements related to CFsubhr files
        for elem in tree.iter(tag='file'):
            if elem.attrib.get('id') is not None:
                if ("_CFsubhr_" in elem.attrib.get('id')):
                    write=True
                    elem.set('output_freq', '2ts')
                    elem.set('split_start_offset', '2ts')
                    if args.verbosity >= 2: print elem.tag, elem.attrib
                    for childElem in elem.iter(tag='field'):
                        if childElem.attrib.get('field_ref') is not None:
                            if childElem.attrib.get('field_ref') not in excluded_field:
                                childElem.set('freq_op', '2ts')
                                if args.verbosity >= 2: print childElem.tag, childElem.attrib
        # Changes related to cth16
        for elem in tree.iter(tag='grid'):
            if elem.attrib.get('id') is not None:
                if elem.attrib.get('id') == 'grid4Dmisr_greordered':
                    write=True
                    if args.verbosity >= 2: print elem.tag, elem.attrib
                    for childElem in elem.iter(tag='axis'):
                        if args.verbosity >= 2: print childElem.tag, childElem.attrib
                        if childElem.attrib.get('axis_ref') is not None:
                            if childElem.attrib.get('axis_ref') == 'cth16':
                                childElem.set('axis_ref', 'cth')
        # Ready
        if write:
            # Write out the results
            if args.verbosity >= 1: print '\nWriting dr2xml_def_xml=','modified.' + inputFile
            try:
                tree.write ('modified.' + inputFile)
            except:
                print "Write error. Please fix so that it can be parsed."
                traceback.print_exc(file=sys.stdout)
                return

if __name__ == '__main__':

    try:
        # Create the top-level parser
        parser = argparse.ArgumentParser(description='XIOS2 xml files tooling and ironsmith')
        subparsers = parser.add_subparsers(description='Dump, check or modify xios xml files')

        # create the parser for the "dump" command
        parser_dump = subparsers.add_parser('dump',help='Dump the xml content without all the xml\'s ironsmith')
        parser_dump.add_argument('-d', '--depth', type=int, default=None, help='How deep do we go. Full tree by default')
        parser_dump.add_argument('file', nargs='+', help='XIOS xml file(s) to dump')
        parser_dump.set_defaults(func=dump)

        # create the parser for the "tsquery" command
        parser_check = subparsers.add_parser('tsquery', help='query timeseries related parameters from an xml file')
        parser_check.add_argument('--file', nargs=1, required=True, help='XIOS xml timeseries_def type')
        parser_check.set_defaults(func=tsquery)
        
        # create the parser for the "check" command
        parser_check = subparsers.add_parser('check', help='Check consistency between field_def and file_def files')
        parser_check.add_argument('--field', nargs=1, required=True, help='XIOS xml field_def type')
        parser_check.add_argument('--file', nargs='+', required=True, help='XIOS xml file_def type')
        parser_check.set_defaults(func=check, correction=False)

        # create the parser for the "modify" command
        parser_check = subparsers.add_parser('modify', help='Will make sure field_def is a superset of file_def')
        parser_check.add_argument('--field', nargs=1, required=True, help='XIOS xml field_def type')
        parser_check.add_argument('--file', nargs='+', required=True, help='XIOS xml file_def type')
        parser_check.set_defaults(func=check, correction=True)
        
        # create the parser for the "modifydr2xmlIPSL" command
        parser_check = subparsers.add_parser('modifydr2xmlIPSL', help='Adapt dr2xml_file to produce output files according to IPSL requirements')
        parser_check.add_argument('--newPath', nargs=1, required=True, help='Directory output for file produced from dr2xml_file type')
        parser_check.add_argument('--JobName', nargs=1, required=False, help='Name of the simulation to be added in the name of file produced from dr2xml_file type')
        parser_check.add_argument('--ModelName', nargs=1, required=False, help='Name of the model to be added in the name of file produced from dr2xml_file type')
        parser_check.add_argument('--ModelVersion', nargs=1, required=False, help='Model version to be added in global attributes of files produced from dr2xml_file type')
        parser_check.add_argument('--ExperimentName',nargs=1, required=False, help='Name of the experiment')
        parser_check.add_argument('--file', nargs='+', required=True, help='XIOS xml dr2xml_file type')
        parser_check.add_argument('--splitLastDate', nargs=1, required=False, help='Date to be used as split_last_date')
        parser_check.set_defaults(func=modifydr2xmlIPSL)

        # create the parser for the "modifyPath" command
        parser_check = subparsers.add_parser('modifyPath', help='Will make sure dr2xml write files where it should be')
        parser_check.add_argument('--newPath', nargs=1, required=True, help='Directory output for file produced from dr2xml_file type')
        parser_check.add_argument('--filePrefix', nargs=1, required=False, help='Prefix to be added to file produced from dr2xml_file type')
        parser_check.add_argument('--file', nargs='+', required=True, help='XIOS xml dr2xml_file type')
        parser_check.set_defaults(func=modifyPath)

        # create the parser for the "modifyDR2XML" command
        parser_check = subparsers.add_parser('modifyDR2XML', help='Fine tune dr2xml files for CMIP6 production')
        parser_check.add_argument('--file', nargs='+', required=True, help='XIOS xml dr2xml_file type')
        parser_check.set_defaults(func=modifyDR2XML)
        
        # create the parser for the "showtime" command
        parser_check = subparsers.add_parser('showtime', help='Just want to make sure you feel good today')
        parser_check.set_defaults(func=showtime)
        
        # Each possible option
        parser.add_argument('-v', '--verbosity', action='count', default=0)

        # Parse the args.
        args = parser.parse_args()

        # And call whatever function was selected
        args.func(args)
    except KeyboardInterrupt:
        print "Shutdown requested...exiting"
    except Exception:
        traceback.print_exc(file=sys.stdout)
    sys.exit(0)
