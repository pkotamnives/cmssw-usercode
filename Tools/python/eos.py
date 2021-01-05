import sys, subprocess, os, fnmatch, re

url = 'root://eosuser.cern.ch/'
global_url = 'root://cms-xrd-global.cern.ch/'

def _popen(cmd, shell=False):
    if type(cmd) == str:
        cmd = cmd.split()
    return subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=shell)

def _system(cmd):
    if type(cmd) == str:
        cmd = cmd.split()
    #print cmd
    return subprocess.call(cmd, stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT) == 0

def fusemount(fn):
    return '/eos/uscms' + storeonly(fn)

def storeonly(fn):
    if fn.startswith('/store'):
        return fn
    if '/store/' not in fn:
        raise ValueError('fn passed %r does not contain "/store/"' % fn)
    return '/store/' + fn.split('/store/')[1]
    
def canon(fn):
    if fn.startswith('root:'):
        return fn
    elif fn.startswith('/store'):
        return url + fn
    else:
        return os.path.realpath(fn)

def quota():
    user = os.environ.get('USER')
    if not user:
        raise ValueError('no env var $USER')
    x = _popen('eos root://cmseos.fnal.gov quota /eos/uscms/store/user/').communicate()[0].split('\n')
    for i, line in enumerate(x):
        if user in line:
            return x[i-1] + '\n' + line
    return "eos quota command didn't work or didn't find user %s" % user

def exists(fn):
    #fn = storeonly(fn)
    return _system('eos %s ls %s' % (url, fn))

def isdir(fn):
    fn = storeonly(fn)
    return _system('eos %s stat -d %s' % (url, fn))

def isfile(fn):
    fn = storeonly(fn)
    return _system('eos %s stat -f %s' % (url, fn))

def size(fn, _re=re.compile(r'Size: (\d+)\s+')):
    fn = storeonly(fn)
    x = _popen('eos %s stat %s' % (url, fn)).communicate()[0].strip()
    mo = _re.search(x)
    if not mo:
        raise ValueError('could not parse eos stat output: %r' % x)
    return int(mo.group(1))

def mkdir(fn):
    fn = storeonly(fn)
    dn = os.path.dirname(fn)
    if exists(dn):
        return True
    return _system('eos %s mkdir -p %s' % (url, dn))

def ls(path):
    fn = storeonly(path)
    fns = _popen(str('eos %s ls %s' % (url, path))).communicate()[0].strip().split()
    return [os.path.join(path, fn) for fn in fns]

def glob(path, pattern):
    fns = ls(path)
    return [fn for fn in fns if fnmatch.fnmatch(os.path.basename(fn), pattern)]

def cp(src, dst):
    src = canon(src)
    dst = canon(dst)
    if not src.startswith(url) and not dst.startswith(url):
        raise ValueError('refusing to use xrdcp for non eos copy: %r -> %r' % (src, dst))
    return _system('xrdcp -s %s %s' % (src, dst))

def rm(fn):
    fn = storeonly(fn)
    return _system('eos %s rm %s' % (url, fn))

def md5sum(fn):
    fn = storeonly(fn)
    cmd = 'xrdcp -s %s%s -' % (url, fn)
    p  = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p2 = subprocess.Popen(('md5sum',), stdin=p.stdout, stdout=subprocess.PIPE)
    p.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    return p2.communicate()[0].split()[0]

def touch(fn):
    fn = storeonly(fn)
    return _system('eos %s touch %s' % (url, fn))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        cmd = sys.argv[1]

        if cmd == 'md5sum':
            for x in sys.argv[2:]:
                if x.startswith('/store'):
                    print '%s  %s' % (md5sum(x), x)

        elif cmd == 'cp' or cmd == 'mv':
            import general

            dest = sys.argv[-1]
            if not dest.startswith('/store'):
                raise ValueError('dest must start with /store')
            if not exists(dest):
                raise IOError("dest %s doesn't exist" % dest)
            if not dest.endswith('/'):
                dest += '/'
            fns = sys.argv[2:-1]
            if not fns:
                raise ValueError('usage: eos.py mv fn1 [fn2...] /store/.../dest')
            print dest
            for fn in fns:
                dest_fn = os.path.join(dest, os.path.basename(fn))
                print fn, '->', dest_fn
                if not os.path.isfile(fn):
                    raise ValueError("files only, can't recurse")
                if exists(dest_fn):
                    raise IOError("dest exists, refusing to clobber")
                if not cp(fn, dest):
                    raise IOError("problem copying?")
                md5src = general.md5sum(fn)
                md5dst = md5sum(dest_fn)
                print '    ', md5src, '->', md5dst
                if md5src != md5dst:
                    raise ValueError("problem with md5sums: src %s dst %s" % (md5src, md5dst))
                if cmd == 'mv':
                    os.remove(fn)
