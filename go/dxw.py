# -*- coding: utf-8 -*-

from django.shortcuts import render
from django.views.decorators import csrf
from pc import Daoxian

info1 = """
 0 -- 1	65.893
 1 -- 2	103.656
 2 -- 3	96.446
 3 -- 4	82.557
 4 -- 5	已知值
 5 -- 6	71.194
 6 -- 0	102.801
 0 -- 7	63.84
 7 -- 8	70.315
 8 -- 3	47.726
 """

info2 = """
∠0 1 2	231°50′15.1″
∠1 0 6	134°25′21.9″
∠1 0 7	62°25′37.2″
∠1 2 3	251°31′27.3″
∠2 3 4	204°32′33.3″
∠2 3 8	288°18′12.4″
∠3 4 5	246°45′15.2″
∠4 3 8	83°45′53.7″
∠4 5 6	244°09′57.8″
∠5 6 0	215°35′53.8″
∠6 0 7	288°00′32″
∠0 7 8	192°59′25.3″
∠7 8 3	156°15′13.5″
"""

info3 = '''已知点：4 x=97.478(m), y=155.682(m)
已知点：5 x=184.833(m), y=90.147(m)
测边误差：1/2000
测角误差：12"'''


# 接收POST请求数据
def cal_dxw(request):
    ctx = {
        'info1': info1, 'info2': info2, 'info3': info3
    }
    if request.POST:
        ctx['params'] = 'producing'
        try:
            ctx['info1'], ctx['info2'], ctx['info3'] = request.POST['info1'], request.POST['info2'], request.POST['info3']
            dxw = Daoxian(ctx['info1'], ctx['info2'], ctx['info3'])
            dxw.check_info()
            dxw.init_params()
            dxw.init_weight()
            dxw.cal()
            ctx['params'] = '<p>params</p>' + dxw.params_to_html_th()
            ctx['edges'] = '<p>edges observation</p>' + dxw.edges_to_html_th()
            ctx['corners'] = '<p>corners observation</p>' + dxw.corners_to_html_th()
        except:
            ctx['params'] = '<p>error input</p>'
            ctx['edges'] = ''
            ctx['corners'] = ''
    return render(request, "dxw_solve.html", ctx)