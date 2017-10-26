def phipsi(selection):
    r = cmd.get_phipsi(selection)
 
    if r is not None:
        keys = r.keys()
        keys.sort()
 
        cmd.feedback('push')
        cmd.feedback('disable','executive','actions')
        for key in keys:
            cmd.iterate("(%s`%d)" % key, "print ' %-5s " + ("( %4.0f, %4.0f ) " % r[key]) + "'%(resn+'-'+resi+' '+chain+':')")
        cmd.feedback('pop')
cmd.extend("phipsi", phipsi)
