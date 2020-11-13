
import os

from django.utils.encoding import smart_str

def handle_uploaded_file(f):

    filename = f.name
    assert filename.endswith(".sdf")

    temp_file = "temp.sdf"
    with open(temp_file, "wb+") as out_file:
        for chunk in f.chunks():
            out_file.write(chunk)

    pass_out_file = os.path.splitext(filename)[0] +\
        "-PASS-out.sdf"

    cmd = "PASS2019toSDF {} {}".format(temp_file, pass_out_file)
    print ("executing command:", cmd)

    ret = os.system(cmd)

    assert ret == 0

    return smart_str(pass_out_file)



