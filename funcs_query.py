import numpy as np

def query(meta,data,idx):
    #dic = {'PLATE':plate, 'MJD':mjd, 'FIBERID':fiber}
    plate, mjd, fiber = meta[idx][3], meta[idx][4], meta[idx][5]
    objid = '{}-{}-{:04d}'.format(plate,mjd,fiber)
    dic = {'OBJID':objid}
    mask = np.ones(len(data), dtype=bool)
    try:
        for key, value in dic.items():
            mask &= (data[key] == value)
        idx = np.where(mask)[0][0]
    except IndexError:
        try:
            ra,dec=meta[idx][0],meta[idx][1]
            dic = {'RA':ra, 'DEC':dec}
            for key,value in dic.items():
                mask &= (data[key]<=value+0.01) & (data[key]>=value-0.01)
            idxs = np.where(mask)[0]
            if len(idxs)>1:
                print('Found more than one object for given coordinates.')
                idx = idxs
            else:
                idx = idxs[0]
        except IndexError:
            print('Could not find object for given coordinates.')
            idx = None
        except Exception as error:
            err = type(error).__name__
            print('Something else went wrong:', err)
            idx = None
    except Exception as error:
        err = type(error).__name__
        print('Something else went wrong:', err)
        idx = None
    return idx


def pmf(meta,idx):
    plate, mjd, fiber = meta[idx][3], meta[idx][4], meta[idx][5]
    return [plate,mjd,fiber]

def identify(meta, data, idx):
    obj = pmf(meta,idx)
    return query(data,*obj)
    
