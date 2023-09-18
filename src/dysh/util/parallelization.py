from threading import Thread
import concurrent

def SingleThread(func):
    """
    Decorator that multithreads the target function
    with the given parameters. Returns the thread
    created for the function
    """
    def wrapper(*args, **kwargs):
        thread = Thread(target=func, args=args)
        thread.start()
        return thread
    return wrapper

def MultiThread(func):
    def wrapper(*args, **kwargs):
        results = {}
        # We can use a with statement to ensure threads are cleaned up promptly
        with concurrent.futures.ThreadPoolExecutor() as executor:

            futures = {executor.submit(func, [i]): idx for idx,
                        i in enumerate(args[0])}

            tasks = len(futures)
            tenth = round(tasks / 10)
            print(f'Formed pool of {tasks} tasks')

            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                i = futures[future]
                try:
                    # store result
                    data = future.result()
                    # check to see if in array form
                    if len(data) == 1:
                        data = data[0]
                    results[i] = data
                except Exception as exc:
                    print('{} generated an exception: {}'.format(
                        args[0][i], exc))

                if tenth != 0 and idx != 0 and idx % tenth == 0:
                    print('{}% Done'.format((idx // tenth) * 10))

        # sort and put in array
        final = []
        for k, v in sorted(results.items()):
            final.append(v)

        return final
    return wrapper