import sys
from bdfs.OneDimBDF import D
from bdfs.OuterMinBDF import OuterMin
from bdfs.OuterMaxBDF import OuterMax


def weighted_sum(weights, power):
    """
    Returns the of all weights^exponent.
    Useful to get just get the sum W, using exp = 1.
    Or with exponent |s_i| / D_v(k) to estimate the threshold
    """
    result = 0
    for i in weights:
        result += (i ** power)
    return result


def weight_scaling(weights, thr_con, depth_vol):
    weight_adjust = (thr_con ** depth_vol) * sum(weights)
    return [weight * weight_adjust for weight in weights]


def filter_by_index(positions, max_set):
    return [[row[i] for i in max_set] for row in positions]


def parse_args():
    '''
    Parses command line argument for the cli

    Returns: a dictionary with argument-value pairs. If an argument has no value it is set to true.

    '''
    args = sys.argv[1:]
    args_dict = {}

    i = 0
    while i < len(args):
        if args[i].startswith('-'):
            key = args[i][1:]
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                value = args[i + 1]
                i += 1
            else:
                value = True
            args_dict[key] = value
        i += 1

    return args_dict


def parse_bdf(bdf):
    '''
    Parses a bdf in string representation

    Args:
        bdf: representation as min, max terms. If a bdf is not given as binary, it is converted to one.

    Returns: A BDF

    '''
    bdf = bdf.strip()
    if bdf.isdigit():
        return D(int(bdf))
    if len(bdf) >= 8:  # Anything else cannot be a bdf. 3 for min/max 2 for brackets 1 for comma 2 for sub_bdfs
        inner_expr = bdf[4:-1]
        sub_expressions = split_top_level(inner_expr)
        sub_formulae = [parse_bdf(sub_expr) for sub_expr in sub_expressions]
        if bdf.startswith('min(') and bdf.endswith(')'):
            return create_bdf_tree(sub_formulae, OuterMin)
        elif bdf.startswith('max(') and bdf.endswith(')'):
            return create_bdf_tree(sub_formulae, OuterMax)
    raise ValueError(f"Invalid bdf: {bdf}")


def split_top_level(expression):
    """ Splits the expression by commas, considering nested parentheses. """
    sub_expressions = []
    start = 0
    depth = 0

    for i, char in enumerate(expression):
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
        elif char == ',' and depth == 0:
            sub_expressions.append(expression[start:i].strip())
            start = i + 1

    sub_expressions.append(expression[start:].strip())
    return sub_expressions


def create_bdf_tree(sub_bdf, outer_type):
    """ Creates a binary tree bdf representation"""
    if len(sub_bdf) == 1:
        return sub_bdf[0]
    else:
        mid = len(sub_bdf) // 2
        return outer_type(
            create_bdf_tree(sub_bdf[:mid], outer_type),
            create_bdf_tree(sub_bdf[mid:], outer_type)
        )


