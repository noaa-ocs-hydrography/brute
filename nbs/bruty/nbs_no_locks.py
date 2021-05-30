EXCLUSIVE, SHARED, NON_BLOCKING = 1, 1, 1


class Lock:
    def __init__(self, fname, mode='r', *args, **kywds):
        self.fname = fname
        self.mode = mode

    def acquire(self, *args, **kywds):
        return True

    def release(self):
        pass

    @property
    def fh(self):
        # make this to mimic the portalocker implementation
        return self.openfile()

    def openfile(self):
        handle = open(self.fname, self.mode)
        return handle

    def __enter__(self):
        return self.openfile()

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


class AreaLock:
    def __init__(self, *args, **kywds):
        pass

    def acquire(self, *args, **kywds):
        return True

    def release(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


class LockNotAcquired(Exception):
    pass
