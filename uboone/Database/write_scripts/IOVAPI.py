import psycopg2
import time
import sys
from datetime import datetime, timedelta
from dbdig import DbDig
from psycopg2 import extensions
import threading

Version = "$Id: IOVAPI.py,v 1.29 2013/06/12 02:31:50 ivm Exp $"

class   IOVCache:
    def __init__(self, n_to_keep = 1):
        self.Cache = {}     # {folder_name:[iovs]}
        self.NToKeep = n_to_keep
        self.Lock = threading.RLock()
        
    def getSnapshot(self, iov, tag, folder):
        self.Lock.acquire()
        folder_name = folder.Name
        found = None
        ifound = None
        
        folder_cache = self.Cache.get(folder_name,[])
        for i, s in enumerate(folder_cache):
            if s.iovid() == iov and s.Tag == tag:
                found = s
                ifound = i
                break
        if found:        
            # temporarily remove from cache
            folder_cache = folder_cache[:ifound] + folder_cache[ifound+1:]
            #print ' found cache'
        else:
            # if not found, call the folder
            self.Lock.release()
            data = folder.getSnapshotByIOV(iov, tag=tag)
            found = data
            self.Lock.acquire()
        if found:
            folder_cache = [found] + folder_cache
            self.Cache[folder_name] = folder_cache[:self.NToKeep]
        self.Lock.release()
        return found
        
    def clear(self):
        self.Cache = {}

class IOVDB:
    def __init__(self, db = None, connstr = None, namespace = 'public',
            role = None,
            cache = None):
        assert not (db == None and connstr == None)
        self.Connstr = connstr
        self.Namespace = namespace
        self.Role = role
        if db:
            self.DB = db
        else:
            self.reconnect()
        self.DataTypes = {}
        self.Cache = cache
        
    def flushCache(self):
        self.Cache.clear()

    def getCached(self, folder, iovid, tag):
        s = None
        if self.Cache:
            s = self.Cache[(folder, iovid, tag)]
        return s
        
    def putInCache(self, folder, iovid, tag, s):
        if self.Cache:
            self.Cache[(folder, iovid, tag)] = s

    def reconnect(self):
        assert self.Connstr != None
        self.DB = psycopg2.connect(self.Connstr)
        self.DB.set_isolation_level(extensions.ISOLATION_LEVEL_AUTOCOMMIT)
        if self.Role:
            self.DB.cursor().execute("set role %s" % (self.Role,)) 
            c = self._cursor()
            c.execute("show role")
            #print "Role %s set, %x" % (c.fetchone()[0],id(self.DB))           
        self.DB.cursor().execute("set search_path to %s" % (self.Namespace,))
        
    def disconnect(self):
        if self.DB != None:
            self.DB.close()
        self.DB = None

    def useDB(self, db):
        self.DB = db

    def openFolder(self, name, columns = '*'):
        return IOVFolder(self, name, columns)
        
    def createFolder(self, name, columns_and_types, drop_existing=False,
                grants={}):
        # grants =  {'username':'r' or 'username':'rw'}
        colnames = [cn for cn, ct in columns_and_types]
        f = self.openFolder(name, colnames)
        f._createTables(columns_and_types, drop_existing, grants)


    def getFolders (self):

        folders =[]

        dbfolders = DbDig(self.DB).tables(self.Namespace)
        if not dbfolders: return []
        
        #Make sure folder has both tables folder_iovs and folder_data.
        for aval in dbfolders:
            if aval[-5:] == '_iovs':
                fn = aval[:-5]
                if fn + "_data" in dbfolders:
                    folders.append(fn)
                
        return folders


    def probe(self):
        c = self.DB.cursor()
        c.execute("select 1")
        return c.fetchone()[0] == 1
        
    def _cursor(self):
        return self.DB.cursor()


    def _tablesExist(self, *tlist):
        tables = DbDig(self.DB).tables(self.Namespace)
        if not tables:  return False
        tables = [t.lower() for t in tables]
        for t in tlist:
            if not t.lower() in tables: return False
        return True
        
    def _columns(self, table):
        columns = DbDig(self.DB).columns(self.Namespace, table)
        return [(c[0], c[1]) for c in columns]
        
    def _typeName(self, ctype):
        if not self.DataTypes.has_key(ctype):
            c = self._cursor()
            c.execute("select pg_catalog.format_type(%d, null)" % (ctype,))
            self.DataTypes[ctype] = c.fetchone()[0]
        return self.DataTypes[ctype]

class IOVSnapshot:
    def __init__(self, iov, data, columns, interval = None, tag = None):
        # data: list of tuples, 1st element of a tuple is channel id
        # validity interval: tuple (t0, t1) or (t0, None)
        # columns: list of tuple element names and their types. First element
        # of the list is always ("channel","bigint")
        self.Data = data
        self.ValidFrom, self.ValidTo = None, None
        if interval != None:   self.ValidFrom, self.ValidTo = interval
        self.Dict = None
        self.Columns = columns
        self.IOVid = iov
        self.Tag = tag
        #print 'end of constructor'


    # to make it look like a dictionary..
    def __getitem__(self, channel):
        #raise 'call to getitem (%s)' % (channel,)
        self._makeDict()
        return self.Dict[channel]
        
    def items(self):
        self._makeDict()
        return self.Dict.items()

    def keys(self):
        self._makeDict()
        return self.Dict.keys()

    def get(self, channel, default):
        self._makeDict()
        return self.Dict.get(channel, default)
        
    def asDict(self):
        self._makeDict()
        return self.Dict
        
    def asList(self):
        return self.Data

    def setInterval(self, t0, t1):
        self.ValidFrom = t0
        self.ValidTo = t1

    def iov(self):  return (self.ValidFrom, self.ValidTo)

    def columns(self):  return self.Columns

    def _makeDict(self):
        if self.Dict == None:
            data = {}
            for tup in self.Data:
                channel = tup[0]
                data[channel] = tup[1:]
            self.Dict = data
        
    def iovid(self):
        return self.IOVid

class IOVFolder:
    def __init__(self, db, name, columns = '*'):
        self.DB = db
        self.Name = name
        self.TableIOVs = '%s_iovs' % (self.Name,)
        self.TableData = '%s_data' % (self.Name,)
        self.TableTags = '%s_tags' % (self.Name,)
        self.TableTagIOVs = '%s_tag_iovs' % (self.Name,)
        self.Columns = columns
        if self.Columns == '*':
            self.Columns = self._getDataColumns()
        self.Colstr = ','.join(self.Columns)

    def exists(self):
        return self.DB._tablesExist(self.TableData, self.TableIOVs)

    def _getDataColumns(self):
        columns = self.DB._columns(self.TableData)
        return [n for n, t in columns if not (n == 'channel' or n.startswith('__'))]

    def getIOVs(self, t, window = 0, tag = None):
        
        db = self.DB
        c = db._cursor()
        t0 = datetime.fromtimestamp(t)
        t1 = t0 + timedelta(seconds = window)
        
        #raise '%s %s' % (t0, t1)
        if tag:

            c.execute("""select i.iov_id, i.begin_time 
                            from %s i, %s ti
                            where i.begin_time between '%s' and '%s' and
                                i.iov_id = ti.iov_id and
                                ti.tag = '%s'
                            order by i.begin_time""" % 
                            (self.TableIOVs, self.TableTagIOVs, t0, t1, tag))
        else:
            c.execute("""select iov_id, begin_time 
                            from %s i
                            where i.begin_time between '%s' and '%s' 
                                    and active
                            order by i.begin_time""" % 
                            (self.TableIOVs, t0, t1))
            
        lst = c.fetchall()
        if tag:
            c.execute("""select i.iov_id, i.begin_time
                from %s i, %s ti
                where begin_time <= '%s'
                    and i.iov_id = ti.iov_id and ti.tag = '%s'
                order by begin_time desc
                limit 1""" % (self.TableIOVs, self.TableTagIOVs, t0, tag))
        else:
            c.execute("""select iov_id, begin_time
                from %s
                where begin_time <= '%s' and active
                order by begin_time desc
                limit 1""" % (self.TableIOVs, t0))
        tup = c.fetchone()
        if tup and (not lst or lst[0][0] != tup[0]):
            lst = [tup] + lst
        return lst

    def getNextIOV(self, t=None, iovid=None, tag=None):
        assert not (t is None and iovid is None)
        db = self.DB
        c = db._cursor()
        if t == None:
            c.execute("""select begin_time 
                            from %s    
                            where iov_id = %s""" % (self.TableIOVs, iovid))
            tup = c.fetchone()
            if tup == None: return None
            t = tup[0]
        #print ' next iov, t= %s' %t
        if tag:
            #print ' next iov, using tag = %s ' %tag
            tagsql ="""select i.iov_id, i.begin_time
                from %s i, %s ti
                where begin_time > '%s' and
                    ti.tag = '%s' and
                    ti.iov_id = i.iov_id
                order by begin_time limit 1""" % (self.TableIOVs, self.TableTagIOVs, t, tag)
            #print ' next iov, sql = \n %s \n' %tagsql

            c.execute("""select i.iov_id, i.begin_time
                from %s i, %s ti
                where begin_time > '%s' and
                    ti.tag = '%s' and
                    ti.iov_id = i.iov_id
                order by begin_time limit 1""" % (self.TableIOVs, self.TableTagIOVs, t, tag))
        else:
            sql = """select iov_id, begin_time
                from %s
                where active and begin_time > '%s'
                order by begin_time
                limit 1""" % (self.TableIOVs, t)
            #print ' next iov, no tag sql= \n %s \n ' %sql
            
            c.execute("""select iov_id, begin_time
                from %s
                where active and begin_time > '%s'
                order by begin_time
                limit 1""" % (self.TableIOVs, t))

        tup = c.fetchone()
        if tup == None: return None, None
        return tup
        
    def getIOV(self, iovid, tag=None):
        c = self.DB._cursor()

        if tag:
            tagsql = """select i.begin_time
                from %s i, %s ti
                where 
                    ti.tag = '%s' and
                    ti.iov_id = i.iov_id and i.iov_id = %s""" % (self.TableIOVs, self.TableTagIOVs,tag,iovid)
            #print ' getIOV, tag sql = %s \n' %tagsql
            c.execute(tagsql)
        
        else:
            sql ="""select begin_time
                from %s
                where iov_id = %s""" % (self.TableIOVs, iovid)
            #print ' getIOV, no tag sql = \n %s \n' %sql
                
            c.execute("""select begin_time
                from %s
                where iov_id = %s""" % (self.TableIOVs, iovid))

        tup = c.fetchone()
        if not tup:
            #print 'getiov , No data found'   
            return None, None
        t0 = tup[0]
        next_iov, t1 = self.getNextIOV(t0, tag=tag)
        #print ' getIOV, t0= %s, t1=%s' %(t0,t1)
        return t0, t1 

    def getIOV_old(self, iovid, tag=None):
        c = self.DB._cursor()
        c.execute("""select begin_time
                from %s
                where iov_id = %s""" % (self.TableIOVs, iovid))
        tup = c.fetchone()
        if not tup:   return None, None
        t0 = tup[0]
        next_iov, t1 = self.getNextIOV(t0, tag=tag)
        #print ' getIOV, t0= %s, t1=%s' %(t0,t1)
        return t0, t1 
        
    def getTuple(self, t, channel=0, tag=None):
        "returns the tuple for given time and single channel"
        data = self.getData(t, tag=tag)
        if data == None:    return None, None
        tup = data.get(channel, None)
        return tup, data.iov()
        
    def getData(self, t = None, iovid = None, chanComps = None, tag=None):
        "returns dictionary {channel -> (data,...)} for given time"
        #print ' in getdata, iovid = %s, chanComps = %s, tag=%s ' %(iovid, str(chanComps), tag)

        assert not (t == None and iovid == None)
        if type(t) == type(1.0) or type(t) == type(1):
            t = datetime.fromtimestamp(t)
        if iovid == None:
            iovid, t0, t1 = self.findIOV(t, tag=tag)
        if iovid == None:
            return None
        
        if chanComps == None:
            s = self.DB.getCached(self.Name, iovid, tag)
            if s == None:
                s = self.getSnapshotByIOV(iovid, tag=tag)
                if s:
                    self.DB.putInCache(self.Name, iovid, tag, s)
        else:
            s = self.getSnapshotByIOV(iovid, chanComps, tag=tag)
        return s or None

    def addData_old(self, t, data, override_future = False):
        if type(t) in (type(1), type(1.0)):
            t = datetime.fromtimestamp(t)
        db = self.DB
        c = db._cursor()
        c.execute("begin")
        if override_future:
            # hope CASCADE will delete data
            c.execute("delete from %s where begin_time >= '%s'" %
                (self.TableIOVs, t))
        iov = self._insertSnapshot(t, data, c)
        #c.execute("commit")

    def addData(self, t, data, combine = False, override_future = 'no'):
        # override future:
        #   'no' or None or '' - do not override
        #   'hide' - mark them inactive
        #   'delete' - delete IOVs
        if type(t) in (type(1), type(1.0)):
            t = datetime.fromtimestamp(t)
        if combine:
            s = self.getData(t)
            old_data = {}
            if s:   old_data = s.asDict()
            old_data.update(data)
            data = old_data
            
        
        db = self.DB
        c = db._cursor()
        
        #c.execute("begin")
        if override_future == 'hide':
            c.execute("update %s set active='false' where active and begin_time >= '%s'" %
                (self.TableIOVs, t))
        elif override_future == 'delete':
            c.execute("""delete from %s i where active and begin_time >= '%s'
                            and not exists (
                                select * from %s ti
                                    where ti.iov_id = i.iov_id)""" %
                (self.TableIOVs, t, self.TableTagIOVs))
        iov = self._insertSnapshot(t, data, c)
        #c.execute("commit")

    def overrideInterval(self, t0, t1, data, mode = 'hide'):
        # mode: 'hide' - mask existing intervals by making them inactive
        #       'delete' - delete existing intervals 
        assert t1 >= t0
        db = self.DB
        c = db._cursor()
        end_data = self.getData(t1)
        self.addData(t0, data)
        if mode == 'delete':
            c.execute("""delete from %s i
                            where   active and
                                    begin_time >= '%s' and
                                    begin_time < '%s' and
                                    not exists (
                                        select * from %s ti
                                        where ti.iov_id = i.iov_id)""" % 
                    (self.TableIOVs, t0, t1, self.TableTagIOVs))
        else:
            c.execute("""update %s
                            set active = 'false'
                            where   active and 
                                    begin_time >= '%s' and
                                    begin_time < '%s'""" % (self.TableIOVs, t0, t1))
        self.addData(t1, end_data)

    def tagIOV(self, iovid, tag):
        c = db._cursor()
        c.execute("""
            insert into %s(tag, iov_id) values ('%s', %s)""" % 
            (self.TableTagIOVs, tag, iovid))
        
    def untagIOV(self, iovid, tag):
        c = db._cursor()
        c.execute("""
            delete 
                from %s
                where tag='%s' and iovid=%s""" % 
            (self.TableTagIOVs, tag, iovid))
        

    def findIOV(self, t, tag=None):
        c = self.DB._cursor()
        if tag:
            c.execute("""
                select i.iov_id, i.begin_time
                    from %s i, %s ti
                    where i.begin_time <= '%s' and
                        ti.iov_id = i.iov_id and
                        ti.tag = '%s'
                    order by i.begin_time desc, i.iov_id desc
                    limit 1""" % (self.TableIOVs, self.TableTagIOVs, t, tag))
        else:
            c.execute("""select iov_id, begin_time
                from %s
                where begin_time <= '%s' and active
                order by begin_time desc, iov_id desc
                limit 1""" % (self.TableIOVs, t))
        x = c.fetchone()
        if not x:
            # not found
            return None, None, None
        iov, t0 = x
        next_iov, t1 = self.getNextIOV(t, tag=tag)
        return iov, t0, t1

    def getSnapshotByIOV(self, iov, chanComps=None, tag=None):
        
        #print ' getSnapshotByIOV iov = %s, tag = %s' %(iov,tag)
        t0, t1 = self.getIOV(iov, tag=tag)
        if t0 == None:  return None
        c = self.DB._cursor()
        
        andSql=''
        if chanComps: print ' chancomps = %s ' %str(chanComps)
        if chanComps:
          for ak in chanComps.keys():
             #print ' snap ,%s, %s' %(ak, chanComps[ak])
             andSql = andSql + " and %s = %s " %(ak,chanComps[ak])
        
        sql = """select channel, %s
                from %s
                where __iov_id = %d %s
                order by channel""" %(self.Colstr, self.TableData, iov, andSql)

        #print 'sql = %s ' %sql
        c.execute(sql)
        desc = c.description
        data = c.fetchall()
        columns = []
        for tup in desc:
            cname, ctype = tup[:2]
            columns.append((cname, self.DB._typeName(ctype)))
        s = IOVSnapshot(iov, data, columns, interval = (t0, t1))
        return s

    def iovsAndTime(self,tagval):
       
        iovData=[]
        iovsTags = []       
        t0=''
        t1=''
        c = self.DB._cursor()
    
        if tagval:
            sql = " select iov_id,tag from %s where tag='%s' "%(self.TableTagIOVs,tagval)
            c.execute(sql)
            tagdata = c.fetchall()
            if tagdata:
                for aval in tagdata:
                    iovid,tag=aval
                    iovsTags.append(iovid)                      
            
        c.execute("select distinct(iov_id),begin_time from %s order by begin_time" %(self.TableIOVs))
        data = c.fetchall()
        if data:
            for aval in data:
                iov = aval[0]
                if tagval and iov not in iovsTags:
                    continue
                    
                t0 =  aval[1]
                next_iov, t1 =  self.getNextIOV(t0)
                iovData.append( (iov,t0,t1))
                
        return iovData
        

    def createTag(self, tag, comments=''):
        c = self.DB._cursor()
        c.execute("""
            insert into %s(tag, created, comments) values('%s', now(), '%s')""" %
            (self.TableTags, tag, comments))
        c.execute("""
            insert into %s(tag, iov_id) 
                (select '%s', iov_id 
                    from %s
                    where active)""" % (self.TableTagIOVs, tag, self.TableIOVs))
        c.execute('commit')
                    

    def tags(self):
        c = self.DB._cursor()
        
        if self.DB._tablesExist(self.TableTags):

            c.execute("""select tag from %s""" % (self.TableTags,))
            tagdata=c.fetchall()
            if tagdata:
                return [x[0] for x in tagdata]
        else:
            #print '\n table %s does not exists..' %self.TableTags
            return None

    def _insertSnapshot(self, t, data, cursor):
        # data is a dictionary {channel -> tuple}
        c = cursor
        #try:
            #sql = "insert into %s(begin_time) values ('%s')" %(self.TableIOVs, t,)

        c.execute("insert into %s(begin_time) values ('%s')" %
                (self.TableIOVs, t,))
        #except:
        #    msg = 'MSG: %s %s<br>SQL statement: %s' % (sys.exc_type, sys.exc_info()[1], sql)
        #    print '%s ' %msg

        c.execute("select lastval()")
        iov = c.fetchone()[0]
        for cn, tup in data.items():
            vals = ['%s' % (iov,), '%s' % (cn,)] + ["'%s'" % (x,) for x in tup]
            vals = ','.join(vals)
            c.execute("insert into %s(__iov_id, channel, %s) values(%s)" %
                (self.TableData, self.Colstr, vals))
        c.execute('commit')
        return iov
        
    def _createTables(self, coltypes, drop_existing, grants={'public':'r'}):
        # grants: {username:'r' or 'rw'}
        db = self.DB
        c = db._cursor()
        columns = ','.join(['%s %s' % (cn, ct) for cn, ct in coltypes])
        #print columns
        if drop_existing:
            try:    
                c.execute("""
                    drop table %(data_table)s;""" %
                    {
                        'tag_iovs_table':self.TableTagIOVs,
                        'tags_table':self.TableTags,
                        'data_table':self.TableData, 
                        'iovs_table':self.TableIOVs,
                    })
            except:
                pass
            try:    
                c.execute("""
                    drop table %(tag_iovs_table)s;""" %
                    {
                        'tag_iovs_table':self.TableTagIOVs,
                        'tags_table':self.TableTags,
                        'data_table':self.TableData, 
                        'iovs_table':self.TableIOVs,
                    })
            except:
                pass
            try:    
                c.execute("""
                    drop table %(tags_table)s;""" % 
                    {
                        'tag_iovs_table':self.TableTagIOVs,
                        'tags_table':self.TableTags,
                        'data_table':self.TableData, 
                        'iovs_table':self.TableIOVs,
                    })
            except:
                pass
            try:    
                c.execute("""
                    drop table %(iovs_table)s;""" % 
                    {
                        'tag_iovs_table':self.TableTagIOVs,
                        'tags_table':self.TableTags,
                        'data_table':self.TableData, 
                        'iovs_table':self.TableIOVs,
                    })
            except:
                pass
        if not self.exists(): 
            #c.execute("show role")
            #print 'Role=%s, %x' % (c.fetchone()[0], id(self.DB.DB))
            sql = """
                create table %(iovs_table)s (
                    iov_id      bigserial   primary key,
                    begin_time  timestamp,
                    active      boolean    default 'true');
                create index %(iovs_table)s_begin_time_inx on %(iovs_table)s(begin_time);

                create table %(tags_table)s (
                    tag         text    primary key,
                    created     timestamp,
                    comments    text
                );

                create table %(tag_iovs_table)s (
                    tag     text    references %(tags_table)s(tag) on delete cascade,
                    iov_id  bigint  references %(iovs_table)s(iov_id) on delete cascade,
                    primary key (tag, iov_id)
                );

                create table %(data_table)s (
                    __iov_id  bigint  references %(iovs_table)s(iov_id) on delete cascade,
                    channel bigint default 0, 
                    %(columns)s,
                    primary key(__iov_id, channel)); """ % \
                    {
                        'data_table':self.TableData, 'iovs_table':self.TableIOVs,
                        'tags_table':self.TableTags, 'tag_iovs_table':self.TableTagIOVs,
                        'columns':columns
                    }
            #print sql
            if grants:
                for u, p in grants.items():
                    sql += """
                        grant select on %s, %s, %s, %s to %s; """ % (
                            self.TableData, self.TableIOVs, self.TableTags, self.TableTagIOVs, u)
                    if p == 'rw':
                        sql += """
                            grant insert, update, delete on %s, %s, %s, %s to %s; """ % (
                                self.TableData, self.TableIOVs, self.TableTags, self.TableTagIOVs, u)
            #print sql
            c.execute(sql)       
                 
if __name__ == '__main__':
    import sys, random
    from datetime import datetime, timedelta
    
    db = IOVDB(connstr=sys.argv[1])
    
    db.createFolder("sample",                       # folder name
        [
            ('x1', 'float'),('x2', 'float'),('x3', 'float'),
            ('x4', 'float'),('x5', 'float'),('x6', 'float'),
            ('x7', 'float'),('x8', 'float'),('x9', 'float'),
            ('n', 'int'),('s', 'text')
        ],  # columns and their types
        drop_existing=True)
    
    f = db.openFolder("sample",     # folder name 
        ['x1','x2','x3','x4','x5','x6','x7','x8','x9','n', 's']              # columns we want to see
        )
    
    channels = 3000
    
    time0 = datetime.now() - timedelta(days=1)  # 1 day ago
    totalt = 0.0
    snapshots = 10
    for s in range(snapshots):               # record 100 meaurements
        t = time0 + timedelta(minutes=s*60)  # every 10 minutes
        t0 = time.time()
        data = {}
        for i in range(channels):
            ic = i * 10
            tup = (random.random(),random.random(),random.random(), 
                random.random(),random.random(),random.random(),
                random.random(),random.random(),random.random(),
                int(random.random()*100), 'channel=%s' % (ic,))
            data[ic] = tup
        f.addData(t, data)
        totalt += time.time() - t0
        
    print 'populated from %s to %s' % (
        time.mktime(time0.timetuple()),time.mktime(t.timetuple()))
        
    print 'Write time per snapshot (sec):', totalt/snapshots
    
    tx = time0 + timedelta(minutes=45)       # t0 + 45 minutes
    t0 = time.time()
    tuple, valid = f.getTuple(tx, 70)
    print "Data for channel=70:", tuple
    print "Valid from %s to %s" % valid
    print "Time:", time.time() - t0
    
    tx = time0 + timedelta(minutes=46)         # t0 + 46 minutes - should get data from cache
    t0 = time.time()
    tuple, valid = f.getTuple(tx, 70)
    print "Data for channel=70:", tuple
    print "Valid from %s to %s" % valid
    print "Time:", time.time() - t0
    
    tx = time0 + timedelta(minutes=25)         # t0 + 25 minutes - should go to the db
    t0 = time.time()
    tuple, valid = f.getTuple(tx, 70)
    print "Data for channel=70:", tuple
    print "Valid from %s to %s" % valid
    print "Time:", time.time() - t0
            
        
        
 
