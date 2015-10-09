class DbDig:

    def __init__(self, conn):
        """
        """
        self.Conn = conn
        self.Dsn = {}
        dsn = self.Conn.dsn
        dsn = dsn.split(' ')
        for v in dsn:
            v = v.split('=')
            if v[0] == 'password':
                continue
            self.Dsn[v[0]] = v[1]

    def dsn(self):
        return self.Dsn

    def dbases(self):
        """
        Find all databases
        """
        c = self.Conn.cursor()
        sql = """select datname from %s""" % ('pg_database',)

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            return [d[0] for d in dd]
        else:
            return None


    def nspaces(self):
        """
        Find all namespaces
        """
        c = self.Conn.cursor()
        sql = """select * from pg_namespace where nspname !~ 'pg_' and nspname != 'information_schema'"""

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            return [d[0] for d in dd]
        else:
            return None


    def tables(self, nspace):
        """
        Find all tables
        """
        c = self.Conn.cursor()
        sql = """select relname from pg_class where relnamespace=(select OID from pg_namespace where nspname='%s') and relkind='r'""" % (nspace,)

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            return [d[0] for d in dd]
        else:
            return None

    def columns(self, nspace, table):
        """
        Find all columns
        """
        c = self.Conn.cursor()
        sql = """SELECT a.attname as "Column",
          pg_catalog.format_type(a.atttypid, a.atttypmod) as "Type",
          (SELECT substring(pg_catalog.pg_get_expr(d.adbin, d.adrelid) for 128)
           FROM pg_catalog.pg_attrdef d
           WHERE d.adrelid = a.attrelid AND d.adnum = a.attnum AND a.atthasdef) as "Modifiers",
          a.attnotnull as "Not NULL", pg_catalog.col_description(a.attrelid, a.attnum) as "Description"
        FROM pg_catalog.pg_attribute a
        WHERE a.attrelid = (select OID from pg_catalog.pg_class where relnamespace=(select OID from pg_namespace where nspname='%s') and relname='%s')
          AND a.attnum > 0 AND NOT a.attisdropped
        ORDER BY a.attnum;
        """ % (nspace, table, )

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            return dd
###            return pp(c, dd)
        else:
            return None


    def indexes(self, nspace, table):
        """
        Find indexes
        """
        c = self.Conn.cursor()
        sql = """SELECT indexname,indexdef FROM pg_indexes WHERE schemaname='%s' AND tablename = '%s'
             """ % (nspace, table, )

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            return dd
        else:
            return None


    def pKey(self, nspace, table):
        """
        Find primary keys
        """
        c = self.Conn.cursor()
        sql = """SELECT constraint_name FROM information_schema.table_constraints WHERE table_schema='%s' AND table_name='%s' AND constraint_type='PRIMARY KEY'
             """ % (nspace, table, )

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            #return dd
             return self.keyDef(nspace, dd[0][0])[-1]    # Return only column list
        else:
            return None


    def fKeys(self, nspace, table):
        """
        Find foreign keys
        """
        c = self.Conn.cursor()
        sql = """SELECT constraint_name FROM information_schema.table_constraints WHERE table_schema='%s' AND table_name='%s' AND constraint_type='FOREIGN KEY'
             """ % (nspace, table, )

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            #return dd
            reply = []
            for kn in dd:
                #reply.append(self.keyDef(nspace, kn[0])[2:])   # return only the first pair - (table_name, column_list) for the specified table itself
                reply.append(self.keyDef(nspace, kn[0]))
            return reply
        else:
            return None



    def referringKeys(self, nspace, table):
        """
        Find all keys
        """
        c = self.Conn.cursor()
        sql = """SELECT constraint_name FROM information_schema.table_constraints WHERE
                 table_schema='%s' AND constraint_type='FOREIGN KEY' AND constraint_name IN 
                (SELECT constraint_name  FROM information_schema.constraint_table_usage WHERE table_schema='%s' AND table_name='%s') ORDER BY table_name
              """ % (nspace, nspace, table, )

        c.execute(sql)
        dd = c.fetchall()
        c.close()
        if dd:
            #return dd
            reply = []
            for kn in dd:
                #reply.append(self.keyDef(nspace, kn[0])[:2])   # Return only the second pair - (table_name, column_list) for the referred table
                reply.append(self.keyDef(nspace, kn[0]))
            return reply
        else:
            return None


    def keyDef(self, nspace, kname):
        """ 
        Find key details
        """
        # Referring table
        c = self.Conn.cursor()
        sql = """SELECT table_name, column_name FROM information_schema.key_column_usage
                WHERE table_schema='%s' AND constraint_name='%s'
              """ % (nspace, kname, )
        c.execute(sql)
        rr = c.fetchall()
        if rr:
          rtname = rr[0][0]
          rcname = tuple([v[1] for v in rr])
        else:
            return None
        # Referred table
        #sql = """SELECT table_name, column_name FROM information_schema.constraint_column_usage
        #        WHERE table_schema='%s' AND constraint_name='%s'
        #      """ % (nspace, kname, )
        sql =   """select table_name, column_name from
                        (SELECT nr.nspname, r.relname, a.attname, c.conname
                          FROM pg_namespace nr, pg_class r, pg_attribute a,
                                pg_namespace nc, pg_constraint c
                                         WHERE nr.oid = r.relnamespace AND r.oid = a.attrelid AND
                                nc.oid = c.connamespace AND
                                               CASE
                                                   WHEN c.contype = 'f'::"char" THEN r.oid =
                                c.confrelid AND (a.attnum = ANY (c.confkey))
                                                   ELSE r.oid = c.conrelid AND (a.attnum = ANY (c.conkey))
                                               END AND NOT a.attisdropped AND (c.contype = ANY
                                (ARRAY['p'::"char", 'u'::"char", 'f'::"char"])) AND r.relkind =
                                'r'::"char") 
                        as x(namespace, table_name, column_name, constraint_name)
                where namespace='%s' and constraint_name='%s'""" % (nspace, kname,)

        c.execute(sql)
        pp = c.fetchall()
        if pp:
          ptname = pp[0][0]
          pcname = tuple([v[1] for v in pp])
        else:
            return None

        c.close()
        # Return combined result
        return (rtname, rcname, ptname, pcname)
        #return (rtname, rcname)


def pp(cursor, data=None, rowlens=0):
    """
    """
    d = cursor.description
    if not d:
        return "#### NO RESULTS ###"
    names = []
    lengths = []
    rules = []
    if not data:
        t = cursor.fetchall()
    for dd in d:    # iterate over description
###     l = dd[1] # Should it be dd[2] i.e. display_size?
        l = dd[2]
        if not l:
            l = 12             # or default arg ...
        l = max(l, len(dd[0])) # handle long names
        names.append(dd[0])
        lengths.append(l)
    for col in range(len(lengths)):
        if rowlens:
            rls = [len(str(row[col])) for row in data if row[col]]
            lengths[col] = max([lengths[col]]+rls)
        rules.append("-"*lengths[col])
    format = " ".join(["%%-%ss" % l for l in lengths])
    result = [format % tuple(names)]
    result.append(format % tuple(rules))
    for row in data:
        result.append(format % row)
    return "\n".join(result)
    #return result


if __name__ == '__main__':
    import psycopg2
    import getopt, sys
    
    opts, args = getopt.getopt(sys.argv[1:], 'p:h:U:W:')

    port     = 5432
    dbname   = 'postgres'  
    user     = 'postgres'
    host     = 'localhost'
    password = ''

    for opt, val in opts:
        if opt=='-p':
            port = val
        if opt=='-U':
            user = val
        if opt=='-h':
            host = val
        if opt=='-W':
            password = val
    if args:
        dbname = args[0]  
    dsn = "dbname='%s' user='%s' host='%s' port=%s password='%s'" % (dbname, user, host, port, password)
###    print "DSN: %s" % "dbname='%s' user='%s' host='%s' port=%s password='************'" % (dbname, user, host, port, )
      
    try:
        conn = psycopg2.connect(dsn)
    except:
        print "I am unable to connect to the database"
        sys.exit(1)


###    c = conn.cursor()

###    c.execute("select * from %s where 1=0" % ('pg_database',))
###    print c.description
###    print
###    c.close()


    dddd = DbDig(conn)

###    print dddd.dsn()

    if len(args)==0:
        dd = dddd.dbases()
        print "Databases:"
        if dd:
            for d in dd:
                print d
        sys.exit(0)


    if len(args)==1:
        dd = dddd.nspaces()
        print "Namespaces:"
        if dd:
            for d in dd:
                print d
        sys.exit(0)


    if len(args)==2:
        nspace = args[1]
        dd = dddd.tables(nspace)
        print "Tables:"
        if dd:
            for d in dd:
                print d
        sys.exit(0)


    if len(args)==3:
        nspace = args[1]
        tname  = args[2]
        dd = dddd.columns(nspace, tname)
        print "Columns:"
        if dd:
            for d in dd:
                print d

        dd = dddd.pKey(nspace, tname)
        print "Primary Key:"
        if dd:
            print dd

        dd = dddd.indexes(nspace, tname)
        print "Indexes:"
        if dd:
            for d in dd:
                print d

        dd = dddd.fKeys(nspace, tname)
        print "Foreign Keys:"
        if dd:
            for d in dd:
                print d

        dd = dddd.referringKeys(nspace, tname)
        print "Referring Keys:"
        if dd:
            for d in dd:
                print d


